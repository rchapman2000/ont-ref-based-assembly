#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2
def helpMessage() {
    log.info"""
ONT Referenced-Based Assembly Pipeline

This pipeline was built specifically for the assembly of Oxford Nanopore Sequencing
reads based on a given reference file.

The pipeline begins by trimming reads using Porechop. The trimmed reads are then aligned to the reference genome using minimap2.
Next, the alignment is corrected using medaka, and variants are called from this correction using Medaka and Longshot.
Finally, depth masking and variants are applied to the reference genome.

USAGE: nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --reference REFERENCE_FASTA --model MEDAKA_MODEL

OPTIONS:

--input INPUT_DIR - [Required] A directory containing demultiplexed, concatenated fastq files generated by Nanopore Sequencing

--output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)

--reference REFERENCE_FASTA - [Required] A reference genome to align reads to.

--model MEDAKA_MODEL - [Required] Medaka requires a 'model' which corresponds to the flow-cell type/basecaller parameters used to correct for errors typical to that technology. Use 'medaka tools list_models' to find this.

OPTIONAL:

    --trim - The pipeline will use Porechop to identify and trim ONT adapters found within reads

    --minReadLen INT - If supplied, the pipeline will perform length filtering using NanoFilt excluding reads less than this size [Default = off]

    --maxReadLen INT - If supplied, the pipeline will perform legnth filtering using NanoFilt excluding reads greater than this size [Default = off]

    --noSecondaryAlignments - If supplied, the pipeline will filter secondary alignments [Default = off]

    --primers PRIMER_BED_FILE - Supply a bed file containing primer coordinates for clipping using 'Samtools Ampliconclip'

    --removeUnclipped - If supplied, the pipeline will remove unclipped reads from a bam alignment (Requires --primers option, [Default = off])

    --minClippedReadLen INT - If primers are supplied via the --primers option, this value will denote the minimum length of a clipped read for it to be considered in downstream analysis (Requires --primers option, [Default = 0])
    
    --minCov INT - The minimum coverage below which a position will be masked [Default = 20]
"""
}

// Function that checks whether a directory ends in a trailing slash.
// This is useful for directory variables that are not parsed into
// file objects in the pipeline (such as the output directory).
def checkDirectoryEnding (fileName) {
    // Grabs the last character in the directory name.
    lastChar = fileName.substring(fileName.length() - 1, fileName.length())

    // Checks whether the last character is slash
    if (lastChar != "/") {
        // If it is not a slash, add that to the directory name.
        fileName = fileName + "/"
    }
    
    // Return the directory name.
    return fileName
}

def createSummaryHeader (trim, minLenFilt, maxLenFilt, primers) {
    FinalHeader = "Sample,Raw Reads,Average Raw Read Length,"

    if (trim != false) {
        FinalHeader = FinalHeader + "Average Trimmed Read Length,"
    }

    if (minLenFilt != 0 || maxLenFilt != 0) {
        FinalHeader = FinalHeader + "Reads Post Length Filter,Average Filtered Read Length,"
    }

    FinalHeader = FinalHeader + "Mapped Reads,"
    
    if (primers != false) {
        FinalHeader = FinalHeader + "Clipped Mapped Reads,"
    }

    FinalHeader = FinalHeader + "SNPs,Indels,Masked Sites,Coverage"

    return FinalHeader
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Defines input parameters. Setting to false by default
// allows us to check that these have been set by the user.
params.input = false
params.reference = false
params.output = false
params.minCov = 20
params.model = false
params.trim = false
params.minReadLen = 0
params.maxReadLen = 0
params.noSecondaryAlignments = false
params.removeUnclipped = false
params.primers = false
params.minClippedReadLen = 0


println "Input Directory: ${params.input}"

include { Setup } from './modules.nf'
include { QC_Report } from './modules.nf'
include { Collect_Raw_Read_Stats } from './modules.nf'
include { Porechop_Trimming } from './modules.nf'
include { QC_Report as QC_Report_Trimmed} from './modules.nf'
include { Length_Filtering } from './modules.nf'
include { QC_Report as QC_Report_Filtered } from './modules.nf'
include { MiniMap2_Alignment } from './modules.nf'
include { Primer_Clipping } from './modules.nf'
include { Medaka_Consensus } from './modules.nf'
include { Call_Variants } from './modules.nf'
include { Generate_Consensus } from './modules.nf'
include { Write_Summary } from './modules.nf'

// Checks the input parameter
if (params.input == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No input directory provided. Pipeline requires an input directory."
    exit(1)
}
else if (!(file(params.input).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.input} is not an existing directory."
    exit(1)
}

// Create a channel for hte input files.
inputFiles_ch = Channel
    // Pull from pairs of files (illumina fastq files denoted by having R1 or R2 in
    // the file name).
    .fromPath("${params.input}*.fastq*")
    // The .fromFilePairs() function spits out a list where the first 
    // item is the base file name, and the second is a list of the files.
    // This command creates a tuple with the base file name and two files.
    .map { it -> [it.getSimpleName(), it]}


// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, ensure that the directory provided ends
    // in a trailing slash (to keep things consistent throughout) the
    // pipeline code.
    outDir = file(params.output).toString()
    println(outDir)
}

// Checks the reference parameter. For this, we cannot use an
// input channel like was used for the input files. Using an input channel
// will cause Nextflow to only iterate once as the reference 
// channel would only only have 1 file in it. Thus, we manually parse
// the reference file into a tuple.
refData = ''
refName = ''
if (params.reference == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: no reference file proivded. Pipeline requires a reference file."
    exit(1)
}
else if (!(file(params.reference).exists())) {
    // If the reference file provided does not exist, notify the user and exit.
    println "ERROR: ${params.reference} does not exist."
    exit(1)
}
else {
    // Process the reference file to be supplied to the index step.
    
    // Parse the file provided into a file object.
    ref = file(params.reference)

    // Grab the basename of the file.
    refName = ref.getBaseName()

    // Place the file basename and file object into
    // a tuple.
    refData = tuple(refName, ref)
}
println(refData)


model = ''
if (params.model == false) {
    println "ERROR: no ONT model provided. Medaka requires a model in the format:"
    println ""
    println "{pore}_{device}_{caller variant}_{caller version}"
    println ""
    println "To see existing models enter: medaka tools list_models"
    exit(1)
}
else {
    model = params.model
}

trimmingEnabled = "DISABLED"
if (params.trim != false) {
    trimmingEnabled = "ENABLED"
}

minReadLenVal = "DISABLED"
minReadLenParam = ""
if (params.minReadLen != 0) {
    minReadLenVal = params.minReadLen + " bp"
    minReadLenParam = "--length ${params.minReadLen}"
}

maxReadLenVal = "DISABLED"
maxReadLenParam = ""
if (params.maxReadLen != 0) {
    maxReadLenVal = params.maxReadLen + " bp"
    maxReadLenParam = "--maxlength ${params.maxReadLen}"
    println maxReadLenParam
}

allowSecondaryAlignVal = "ENABLED"
secondaryAlignParam = "--secondary=yes"
if (args.noSecondaryAlignments != false) {
    secondaryAlignParam = "--secondary=no"
    allowSecondaryAlignVal = "DISABLED"
}

primerFile = ''
primerFileName = 'NONE'
if (params.primers != false) {
    if (file(params.primers).exists()) {
        primerFile = file(params.primers)
        if (primerFile.getExtension() != "bed") {
            println "ERROR: The primer file ${params.primers} is not a .bed file. Please supply a file in BED format."
            exit(1)
        }

        primerFileName = primerFile.getName()
    }
    else {
        println "ERROR: ${params.primers} does not exist. Please provide an existing file."
        exit(1)
    }
}

minClippedReadLenVal = "DISABLED"
if (params.minClippedReadLen != 0 && params.primers == false) {
    println "ERROR: the --minClippedReadLen parameter requires the --primers parameter to be supplied. Please adjust the parameters."
    exit(1)
}
else if (params.minClippedReadLen != 0 && params.primers != false) {
    minClippedReadLenVal = params.minClippedReadLen + " bp"
}

removeUnclippedVal = "DISABLED"
removeUnclippedParam = ""
if (params.removeUnclipped != false && params.primers == false) {
    println "ERROR: the --removeUnclipped parameter requires the --primers parameter to be supplied. Please adjust the parameters."\
    exit(1)
}
else if (params.removeUnclipped != false && params.primers != false) {
    removeUnclippedVal = "ENABLED"
    removeUnclippedParam = "--clipped"
}

summaryHeader = createSummaryHeader(params.trim, params.minReadLen, params.maxReadLen, params.primers)

workflow {

    Setup ( trimmingEnabled, minReadLenVal, maxReadLenVal, refName, allowSecondaryAlignVal, primerFileName, minClippedReadLenVal, removeUnclippedVal, params.minCov, model, summaryHeader, outDir )

    QC_Report( inputFiles_ch, "Raw-Reads", outDir )

    Collect_Raw_Read_Stats( inputFiles_ch )

    if (params.trim != false) {
        Porechop_Trimming( Collect_Raw_Read_Stats.out[0], outDir, Collect_Raw_Read_Stats.out[1] )

        QC_Report_Trimmed( Porechop_Trimming.out[0], "Trimmed-Reads", outDir )

        if (params.minReadLen != 0 || params.maxReadLen != 0) {
            Length_Filtering( Porechop_Trimming.out[0], minReadLenParam, maxReadLenParam, outDir, Porechop_Trimming.out[2] )

            QC_Report_Filtered( Length_Filtering.out[0], "Length-Filtered-Reads", outDir )

            MiniMap2_Alignment( Length_Filtering.out[0], outDir, refData, secondaryAlignParam, Length_Filtering.out[1] )
        }
    }
    else if (params.minReadLen != 0 || params.maxReadLen != 0) {
        Length_Filtering( Collect_Raw_Read_Stats.out[0], minReadLenParam, maxReadLenParam, outDir, Collect_Raw_Read_Stats.out[1] )

        QC_Report_Filtered( Length_Filtering.out[0], "Length-Filtered-Reads", outDir )

        MiniMap2_Alignment( Length_Filtering.out[0], outDir, refData, secondaryAlignParam, Length_Filtering.out[1] )
    }
    else {
        MiniMap2_Alignment( Collect_Raw_Read_Stats.out[0], outDir, refData, secondaryAlignParam, Collect_Raw_Read_Stats.out[1] )
    }

    if (params.primers != false ) {
        Primer_Clipping( MiniMap2_Alignment.out[0], primerFile, removeUnclippedParam, params.minClippedReadLen, baseDir, refData, outDir, MiniMap2_Alignment.out[1] )

        Medaka_Consensus( Primer_Clipping.out[0], model, outDir, Primer_Clipping.out[1] )
    }
    else {
        Medaka_Consensus( MiniMap2_Alignment.out[0], model, outDir, MiniMap2_Alignment.out[1] )
    }

    Call_Variants( Medaka_Consensus.out[0], baseDir, outDir, refData, params.minCov, Medaka_Consensus.out[1] )

    Generate_Consensus( Call_Variants.out[0], baseDir, outDir, refData, params.minCov, Call_Variants.out[2] )

    Write_Summary( Generate_Consensus.out[2], outDir )
}