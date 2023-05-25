// Creates a parameters file and a summary file to 
// be added to later
process Setup {
    input:
        val trimming

        val minReadLen

        val maxReadLen
        // The name of the reference supplied (for use
        // in the parameters file)
        val refName

        val secondaryAlignVal

        val supplementaryAlignVal

        val primerFile

        val minClippedReadLen

        val removeUnclipped
        // The minimum coverage below which a site is masked (for
        // use in the parameters file)
        val minCov
        // The provided medaka model
        val model
        // The header to write to the summary file.
        val summaryHeader
        // The output directory to be used.
        val outDir
        
    output:
        // The parameters file created.
        file "analysis-parameters.txt"
        // The blank summary file to be added to.
        file "stats-summary.csv"

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Creates a parameters file (in case the user needs to look back at how they ran the analysis)
    as well as a blank summary file, which will later be populated with data from each
    sample.

    The parameters file contains:
        1. Whether Porechop trimming was enabled
        2. The minimum read length cutoff
        3. The maximum read length cutoff
        4. The reference file supplied
        5. The primer file supplied (if supplied)
        6. The minimum clipped read length cutoff (if primers are supplied)
        8. Whether unclipped reads are allowed (if primers are supplied)
        9. The medaka model supplied
        10. The minimum base coverage threshold to

    The summary file will always contain:
        1. The sample
        2. Raw Reads
        3. Average Raw Read Length
        4. Average Trimmed Read Length
        5. Mapped reads
        6. SNPs passing filtering
        7. Indels passing filtering
        8. sites masked
        9. coverage
    */
    """
    #!/bin/bash

    touch analysis-parameters.txt

    echo "Porechop Adapter Trimming: ${trimming}" >> analysis-parameters.txt
    echo "Minimum Read Length Cutoff : ${minReadLen}" >> analysis-parameters.txt
    echo "Maxmimum Read Length Cutoff : ${maxReadLen}" >> analysis-parameters.txt
    echo "Reference Supplied : ${refName}" >> analysis-parameters.txt
    echo "Allow Secondary Alignments : ${secondaryAlignVal}" >> analysis-parameters.txt
    echo "Allow Supplementary Alignments: ${supplementaryAlignVal}" >> analysis-parameters.txt
    echo "Primer File Supplied : ${primerFile}" >> analysis-parameters.txt
    echo "Minimum Clipped Read Length : ${minClippedReadLen}" >> analysis-parameters.txt
    echo "Unclipped Reads Removal: ${removeUnclipped}" >> analysis-parameters.txt
    echo "Medaka Model : ${model}" >> analysis-parameters.txt
    echo "Minimum Base Coverage Threshold: ${minCov}" >> analysis-parameters.txt

    touch stats-summary.csv

    echo "${summaryHeader}" > stats-summary.csv
    """
}

process QC_Report {
    input:
        tuple val(base), file(reads)

        val outDir

    output:
        file("Raw-Reads-QC.zip")

    publishDir "${outDir}/${base}-Intermediate-Files/", mode: "copy"

    script:
    """
    #!/bin/bash

    mkdir Raw-Reads-QC

    if [[ -s ${reads} ]]; then
        NanoPlot --fastq ${reads} -o Raw-Reads-QC
    fi

    zip -r Raw-Reads-QC.zip Raw-Reads-QC
    """
}

process QC_Report_Filtered {
    input:
        tuple val(base), file(reads)

        val outDir

    output:
        file("Post-Filtering-QC.zip")

    publishDir "${outDir}/${base}-Intermediate-Files/", mode: "copy"

    script:
    """
    #!/bin/bash

    mkdir Post-Filtering-QC

    if [[ -s ${reads} ]]; then
        NanoPlot --fastq ${reads} -o Post-Filtering-QC
    fi

    zip -r Post-Filtering-QC.zip Post-Filtering-QC
    """
}

process Collect_Raw_Read_Stats {
    input:
        tuple val(base), file(reads)

    output:
        tuple val(base), file(reads)

        env summary

    // The publishDir directive is ignored here because we want the reads
    // to be an output of this process, but to save space, we do not want
    // to make a copy of them.
    script:
    """
    #!/bin/bash

    raw_reads=\$((\$(zcat -f ${reads} | wc -l)/4)) 

    avg_raw_read_len=\$(bioawk -c fastx '{ totalBases += length(\$seq); totalReads++} END{print totalBases/totalReads}' ${reads})

    summary="${base},\$raw_reads,\$avg_raw_read_len"
    """
}

// Detects and trims adapter/barcode sequences from reads using 
// Porechop
process Porechop_Trimming {
    input:
        // Tuple contains the sample basename
        // and the reads to be trimmed
        tuple val(base), file(reads)
        // The name of the output directory
        val outDir

        val existingSummary

    output:
        // Tuple contains the sample basename and the 
        // trimmed reads in fastq format
        tuple val(base), file("${base}-trimmed.fastq.gz")
        // The trimming report generated by porechop
        file "${base}-porechop-report.txt.gz"
        // The summary string
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/Processed-Reads", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${outDir}/${base}-Intermediate-Files/Processed-Reads", mode: 'copy', pattern: "*.txt.gz"

    /*
    Uses Porechop to detect and trim any adapters/barocode
    sequences present within the reads.

    As well, it calculates metrics such as the number of raw
    and trimmed reads as well as the average length of the raw
    and trimmed reads.
    */
    script:
    """
    #!/bin/bash

    porechop -i ${reads} -o ${base}-trimmed.fastq --verbosity 1 > ${base}-porechop-report.txt

    reads_post_trimming=\$((\$(zcat -f ${base}-trimmed.fastq | wc -l)/4))

    avg_trimmed_read_len=\$(bioawk -c fastx '{ totalBases += length(\$seq); totalReads++} END{print totalBases/totalReads}' ${base}-trimmed.fastq)

    gzip ${base}-trimmed.fastq
    gzip ${base}-porechop-report.txt

    summary="${existingSummary},\$reads_post_trimming,\$avg_trimmed_read_len"
    """
}

process Length_Filtering {
    input:
        tuple val(base), file(reads)

        val minLenParam

        val maxLenParam

        val outDir

        val existingSummary

    output:
        tuple val(base), file("${base}-length-filtered.fastq.gz")

        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/Processed-Reads", mode: 'copy', pattern: "*.fastq.gz"

    script:
    """
    #!/bin/bash

    zcat -f ${reads} | chopper ${minLenParam} ${maxLenParam} > ${base}-length-filtered.fastq

    num_filt_reads=\$((\$(cat ${base}-length-filtered.fastq | wc -l)/4)) 

    avg_filt_read_len=\$(bioawk -c fastx '{ totalBases += length(\$seq); totalReads++} END{print totalBases/totalReads}' ${base}-length-filtered.fastq)

    gzip ${base}-length-filtered.fastq

    summary="${existingSummary},\$num_filt_reads,\$avg_filt_read_len"
    """

    //NanoFilt ${minLenParam} ${maxLenParam} ${reads} > ${base}-length-filtered.fastq
}

// Generates an alignment of the asembly contigs to a reference genome
// using minimap.
process MiniMap2_Alignment {
    input:
        // Tuple contains the sample basename
        // and the reads to be aligned.
        tuple val(base), file(reads)
        // The name of the output directory
        val outDir
        // Tuple contains the reference basename 
        // and the reference fasta file to be aligned to.
        tuple val(refName), file(ref)

        val secondaryAlignParam

        val supplementaryAlignParam
        // The existing summary string to be added to.
        val existingSummary

    output:
        // Tuple contains the file basename and the alignment bam file.
        tuple val(base), file("${base}-align.bam")
        // The summary string containing sample name, raw reads, and mapped reads
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/", mode: 'copy', pattern: "${base}-align.bam"
    
    script:
    /*
    Uses minimap2 to align the contigs to the reference fasta.

    Then uses samtools to conver the alignment sam into bam format
    (samtools view). The bam format is then sorted and stored in a bam
    file (samtools sort).
    */
    """
    #!/bin/bash

    minimap2 -ax map-ont ${secondaryAlignParam} ${ref} ${reads} > align.sam

    samtools view -b ${supplementaryAlignParam} align.sam | samtools sort > ${base}-align.bam

    mapped_reads=\$(samtools view -F 0x04 -c ${base}-align.bam)

    average_read_depth=\$(samtools depth -a -J -q 0 -Q 0 ${base}-align.bam | awk -F'\t' 'BEGIN{totalCov=0} {totalCov+=\$3} END{print totalCov/NR}')

    summary="${existingSummary},\$mapped_reads,\$average_read_depth"
    """
}

process Samtools_Primer_Clipping {
    input:
        tuple val(base), file(bam)

        file primers

        val removeUnclippedParam

        val minLen

        val baseDir

        // Tuple contains the reference basename 
        // and the reference fasta file to be aligned to.
        tuple val(refName), file(ref) 

        val outDir

        val existingSummary

    output:
        tuple val(base), file("${base}-clipped-sorted.bam")

        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/", mode: 'copy', pattern: "${base}-clipped-sorted.bam"

    script:
    """
    #!/bin/bash

    samtools ampliconclip -b ${primers} ${removeUnclippedParam} --both-ends --filter-len ${minLen} ${bam} | samtools sort > ${base}-initial-clipped.bam

    python3 ${baseDir}/scripts/find_primer_coverage_region.py -i ${primers} -o ${base}

    bioawk -c fastx '{print \$name"\t0\t"length(\$seq)}' ${ref} > ${ref}-all-sites.bed

    bedtools subtract -a ${ref}-all-sites.bed -b ${base}-coverage-region.bed > ${ref}-uncaptured-regions.bed

    samtools ampliconclip -b ${ref}-uncaptured-regions.bed --both-ends ${base}-initial-clipped.bam | samtools sort > ${base}-clipped-sorted.bam

    clipped_mapped_reads=\$(samtools view -F 0x04 -c ${base}-clipped-sorted.bam)
    
    average_read_depth_post_clipping=\$(samtools depth -a -J -q 0 -Q 0 ${base}-clipped-sorted.bam | awk -F'\t' 'BEGIN{totalCov=0} {totalCov+=\$3} END{print totalCov/NR}')

    summary="${existingSummary},\$clipped_mapped_reads,\$average_read_depth_post_clipping"
    """
}

process XGen_Primer_Clipping {
    input:
        tuple val(base), file(bam)

        file primers

        val outDir

        val existingSummary
    
    output:

        tuple val(base), file("${base}-clipped-sorted.bam")

        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/", mode: 'copy', pattern: "${base}-clipped-sorted.bam"

    script:
    """
    #!/bin/bash

    samtools view -bS ${bam} | samtools sort -n -O sam > ${base}-sort.sam

    primerclip -s ${primers} ${base}-sort.sam ${base}-clip.sam

    samtools view -b ${base}-clip.sam | samtools sort > ${base}-clipped-sorted.bam

    clipped_mapped_reads=\$(samtools view -F 0x04 -c ${base}-clipped-sorted.bam)

    average_read_depth_post_clipping=\$(samtools depth -a -J -q 0 -Q 0 ${base}-clipped-sorted.bam | awk -F'\t' 'BEGIN{totalCov=0} {totalCov+=\$3} END{print totalCov/NR}')

    summary="${existingSummary},\$clipped_mapped_reads,\$average_read_depth_post_clipping"
    """

}

process Medaka_Consensus {
    maxForks 1 

    input:
        // Tuple contains the sample base name and the clipped bam file.
        tuple val(base), file(bam)
        // The medaka model to use (supplied by the user)
        val model
        // The output directory
        val outDir
        // The existing summary string to be added to.
        val existingSummary

    output:
        // Tuple contains the sample name, the corrected consensus by medaka,
        // and the clipped bam file.
        tuple val(base), file("${base}.hdf"), file(bam)
        // The summary string.
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/", mode:'copy', pattern: "${base}.hdf"

    script:
    """
    #!/bin/bash

    samtools index ${bam}

    medaka consensus ${bam} --model ${model} ${base}.hdf --debug

    summary="${existingSummary}"
    """
}

process Call_Variants {
    input:
        // Tuple contains the file basename and alignment bam file
        tuple val(base), file(hdf), file(bam)
        // The script base directory name (to call python scripts)
        val baseDir
        // The output directory name
        val outDir
        // Tuple contains the reference name and reference fasta file
        tuple val(refName), file(ref)
        // The minimum coverage cutoff
        val minCov
        // The existing summary string.
        val existingSummary

    output:
        // Tuple contains the file basename, alignment bamfile, filtered snp vcf, and filtered indel vcf
        tuple val(base), file(bam), file("${base}-snps-filtered.vcf"), file("${base}-indels-filtered.vcf")
        // The VCF file produced by longshot
        file "${base}-longshot.vcf"
        // The summary string with the number of snps and indels added.
        env summary

    publishDir "${outDir}/${base}-Intermediate-Files/", mode: 'copy', pattern: "*.vcf"

    // Need to fix python3 ${baseDir}/scripts/fix_multi_allelic.py -i ${base}-longshot.vcf -o ${base}-biallelic.vcf
    script:
    """
    #!/bin/bash

    samtools index ${bam}

    medaka variant ${ref} ${hdf} ${base}-medaka.vcf

    bgzip ${base}-medaka.vcf
    tabix ${base}-medaka.vcf.gz

    longshot -P 0 -A -a 0 --min_mapq 0 --no_haps -v ${base}-medaka.vcf.gz --bam ${bam} --ref ${ref} --out ${base}-longshot.vcf

    vcftools --keep-only-indels --vcf ${base}-longshot.vcf --recode --recode-INFO-all --stdout > ${base}-indels.vcf
    bgzip ${base}-indels.vcf
    tabix ${base}-indels.vcf.gz

    vcftools --remove-indels --vcf ${base}-longshot.vcf --recode --recode-INFO-all --stdout > ${base}-snps.vcf
    bgzip ${base}-snps.vcf
    tabix ${base}-snps.vcf.gz

    bcftools view -i "((INFO/AC[0] + INFO/AC[1]) >= ${minCov}) && ((INFO/AC[1] / INFO/DP) > (INFO/AC[0] / INFO/DP))" ${base}-indels.vcf.gz > ${base}-indels-filtered.vcf
    num_indels=\$(grep -v "^#" ${base}-indels-filtered.vcf | wc -l)

    bcftools view -i "((INFO/AC[0] + INFO/AC[1]) >= ${minCov}) && ((INFO/AC[1] / INFO/DP) > (INFO/AC[0] / INFO/DP))" ${base}-snps.vcf.gz > ${base}-snps-filtered.vcf
    num_snps=\$(grep -v "^#" ${base}-snps-filtered.vcf | wc -l)

    summary="${existingSummary},\$num_snps,\$num_indels"
    """
}

process Generate_Consensus {
    input:
        // Tuple contains the file basename, the alignment bam, the snp vcf file
        // and the indel vcf file
        tuple val(base), file(bam), file(snps), file(indels)
        // The name of the base directory
        val baseDir
        // The name of the base directory
        val outDir
        // Tuple contains the reference file name and reference file
        tuple val(refName), file(ref)
        // The minimum coverage threshold
        val minCov
        // The existing summary string.
        val existingSummary

    output:
        // Tuple contains the consensus fasta and the sites that were masked in 
        // a bed file.
        tuple val(base), file("${base}-consensus.fasta")
        // The bed file containing the masked sites.
        file "${base}-mask-sites.bed"
        // The summary string with the number of masked positions and coverage
        // added.
        env summary

    publishDir "${outDir}", mode: 'copy', pattern: "${base}-consensus.fasta"
    publishDir "${outDir}/${base}-Intermediate-Files", mode: 'copy', pattern: "${base}-mask-sites.bed"

    script:
    /*
    The script first computes sites to mask by identifying sites that have less than
    the minimum coverage provided. However, there are formatting issues with the pileup
    format that make this difficult. Samtools mpileup's output has the depth in 0-based
    format, which makes it impossible to distinguish between a site with 0 and 1 coverage.

    Thus, the pipeline instead makes use of bedtools subtract. It first creates a pileup for only sites with
    coverage, uses an in-house script to filter sites with less than the minimum coverage
    and converts these into a bed file.

    Next, the pipeline creates a pileup containing every site, and converts this into a bed file using
    the in-house script. 

    Finally, the sites we want to keep (those above the minimum coverage threshold) are substracted
    from the bed file with every site, to give us the low-coverage sites.

    Additionally, there is an interesting case when the sites that fall within a deletion are marked as masked.
    Because masking is applied first, this will cause an error when applying the variants (as the deletion site will contain 
    an N character and will not match the VCF reference). Thus, the pipeline uses bcftools to create a bed file for all indel sites,
    and then subtracts these from the low coverage sites. Now, this results in the sites to mask.

    The bedtools maskfasta command is then used to mask the reference at these positions.

    Then the variants are applied to the mask fasta. The reason this is done after masking is 
    because the pileup (and therefore masking) positions do not account for indels, which would
    shift the genomic coordinates (we would end up masking things we did not want to).

    Finally, the fasta is wrapped to make it visually appealing. 

    Old Code:
    samtools mpileup --no-BAQ -d 100000 -x -A -a -Q 0 -f ${ref} ${bam} > all-sites.pileup
    python3 ${baseDir}/scripts/pileup_to_bed.py -i all-sites.pileup -o all-sites.bed 
    */
    """
    #!/bin/bash

    samtools mpileup --no-BAQ -d 100000 -x -A -Q 0 -f ${ref} ${bam} > ${base}.pileup
    python3 ${baseDir}/scripts/pileup_to_bed.py -i ${base}.pileup -o passed-sites.bed --minCov ${minCov}

    bioawk -c fastx '{print \$name"\t0\t"length(\$seq)}' ${ref} > all-sites.bed

    if [[ -s all-sites.bed ]]; then

        bedtools subtract -a all-sites.bed -b passed-sites.bed > ${base}-low-cov-sites.bed

        bcftools query -f'%CHROM\t%POS0\t%END\n' ${indels} > indel-sites.bed

        bedtools subtract -a ${base}-low-cov-sites.bed -b indel-sites.bed > ${base}-mask-sites.bed

    else
        bioawk -c fastx '{print \$name"\t0\t"length(\$seq)}' ${ref} > ${base}-mask-sites.bed
    fi

    num_mask=\$(bioawk -c bed 'BEGIN{SITES=0} {SITES+=\$end-\$start } END{print SITES}' ${base}-mask-sites.bed)

    bedtools maskfasta -fi ${ref} -bed ${base}-mask-sites.bed -fo masked.fasta

    bgzip ${snps}
    tabix ${snps}.gz

    bgzip ${indels}
    tabix ${indels}.gz

    bcftools consensus -f masked.fasta ${snps}.gz > with-snps.fasta

    bcftools consensus -f with-snps.fasta ${indels}.gz > with-indels-snps.fasta

    bioawk -c fastx '{ gsub(/\\n/,"",seq); print ">${base}-"\$name; print \$seq }' with-indels-snps.fasta > ${base}-consensus.fasta

    seq_len=\$(bioawk -c fastx 'BEGIN{bases=0} { bases+=length(\$seq) } END{print bases}' < ${base}-consensus.fasta)

    coverage=\$(python3 ${baseDir}/scripts/calculate_genome_coverage.py -i ${base}-consensus.fasta)

    summary="${existingSummary},\$num_mask,\$coverage"
    """
}


// Writes a line to the summary file for the sample.
process Write_Summary {
    input:
        // Tuple contains the sample basename and forward/reverse reads (the basename
        // is the only value important to this function).
        val summary
        // The output directory.
        val outDir

    script:
    /*
    The summary string containing the statistics collected as the pipeline
    was run are appended to the summary file.
    */
    """
    #!/bin/bash

    echo "${summary}" >> ${outDir}/stats-summary.csv
    """  

}