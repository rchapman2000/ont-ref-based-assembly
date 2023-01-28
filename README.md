# Nextflow ONT Reference Based Assembly Pipeline

This pipeline automates the process of performing a reference based assembly on Oxford Nanopore (ONT) sequencing data. It also includes parameters to support the processing of amplicon sequencing data on ONT.

## Technical Considerations

### Medaka model

Medaka requires information about the pore, sequencing device, and basecaller. This information is specified to the tool through a 'model', which is a string of text in one of the following formats:

For flowcell versions < 10.4
```
{pore}_{device}_{basecaller accuracy}_{basecaller version} #Example: r941_min_hac_g507
```
or for flowcell versions >= 10.4
```
{pore}_{motor protein version}_{pore speed}_{basecaller accuracy}_{basecaller version} #Example: r1041_e82_400bps_hac_g615
```

The basecaller setting can contain one of three values based on the accuracy setting of the base caller:
- fast - Fast basecalling
- hac - High-accuracy basecalling
- sup - Super-accurate basecalling

As well, we have seen better variant calling accuracy by adding the **_variant_** value to the medaka model between the basecaller accuracy and basecaller version values (Example: r941_min_hac_variant_g507)

Finally, the basecaller version should be the highest listed by medaka that is **less than or equal to** the version of guppy used to base-call the data (Example: if guppy 6.0 was used, g507 should be used).

To see the complete list of medaka models used, enter the command ```medaka tools list_models```


## Installation
To install this pipeline, enter the following commands:
```
# Clone the repository
git clone https://github.com/rchapman2000/ont-ref-based-assembly.git

# Create a conda environment using the provided environment.yml file
conda env create -f environment.yml

# Activate the conda environment
conda activate ONT-RefBasedAssembly
```
## Updating the Pipeline
If you already have the pipeline installed, you can update it using the following commands:
```
# Navigate to your installation directory
cd ont-ref-based-assembly

# Use git to pull the latest update
git pull

# Activate the conda environment and use the environment.yml file to download updates
conda activate ONT-RefBasedAssembly
conda env update --file environment.yml --prune
```

## Usage
To run the pipeline, use the following command:
```
# You must either be in the same directory as the main.nf file or reference the file location.
nextflow run main.nf [OPTIONS] --input INPUT_DIR --output OUTPUT_DIR --reference --model MEDAKA_MODEL
```

### Optional Arguments
The pipeline also supports the following optional arguments:

| Option | Type | Description |
|---|---|---|
| --trim | *None* | If supplied, the pipeline will trim ONT adapters/barcodes from reads using Porechop. [Default = off] |
| --minReadLen | *Int* | If supplied, the pipeline will perform length filtering using Nanofilt excluding reads shorter than this value. [Default = off] |
| --maxReadLen | *Int* | If supplied, the pipeline will perform length filtering using Nanofilt excluding reads larger than this value. [Default = off] |
| --noSecondaryAlignments | *None* | If supplied, the pipeline will only consider primary read alignments. [Default = off] |
| --primers | *File* | The pipeline will take a bed file containing coordinates and perform primer clipping using Samtools. the start and end regions of reads that fall within these coordinates will be soft-clipped from the alignmetn [Default = disabled] |
| --removeUnclipped | *None* | If supplied (and a primer file supplied), the pipeline will remove unclipped reads from the alignment. (Requires --primers option, [Default = off]) |
| --minClippedReadLen | *Int* | If supplied (and a primer file supplied), the pipeline will remove clipped reads shorter than this length from the alignment. (Requires --primers option, [Default = off]) |
| --minCov | *Int* | The minimum depth of coverage required for a position to be called. Sites with a depth less than this value will be reported as "N" [Default = 20] |

To view the list of options from the command line, use the following command:
```
nextflow run main.nf --help
```
