# MHConstructor

## Project Description

Description here

## Installation Guide
To obtain the project directory, the repository can be cloned using the following command:

`git clone https://github.com/wadekj/MHConstructor.git`

To ensure a reliable run, we have containerized our image in [Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/installation.html). We have tested it on Singularity version 3.11.4. The Singularity recipe file, `MHConstructor/container/mhconstructor.def`, can be built using the following command:

`sudo singularity build mhconstructor.sif mhconstructor.def`

This will create a .sif file that will need to be included when you're running the pipeline through Singularity.


## Usage Guide
### To run our pipeline from start to finish with default settings, you can run the following line:

`singularity exec --bind <fastq_location> container/mhconstructor.sif /bin/bash MHCgenerate_v5.sh testID.txt`

### To run specific submodules of interest, you can run the following(s):
- Generate BMH

    `$ python assignMHChaps.py <assignMHChaps.py> 
<HLADRBgenotypes.csv> <C4genotypes.csv>`
- MHC read QC

    `$ sh readQC.sh <sampleID> <path/to/R1> <path/to/R2> 
<path/to/assembly/dir> <#ofThreads>`    
- MHC <i>de novo </i>assembly

    `$ sh run_athena_v3.sh <sampleID>  <path/to/R1> <path/to/R2> 
<path/to/assembly/dir> <BMHap1> <path/to/assembly/dir> 
<path/to/bin/dir> <path/to/samtools> <path/to/Picard> 
<startingReadCount> <insertSize> <kmerSize> 
<expectedAvgCoverage> <BMHap2> <path/to/bin/dir> 
<path/to/Bowtie2>`
- MHC <i>de novo </i>scaffolding

    `$ sh 
refGuidedDeNovoAssembly_velvet_athena_scaffold_v3_target_haploid.s
h <sampleID> <path/to/R1> <path/to/R2> <path/to/assembly/dir> 
<path/to/assembly/dir> <insertSize> <path/to/bin/dir> <BMHap1> 
<BMHap2> <path/to/Bowtie2> <path/to/samtools> <path/to/Picard> 
<startingReadCount>`
- MHC scaffold orientation

    `$ sh orderAssembly.sh <path/to/scaffold/sequence.fasta> 
<BMH1orBMH2> <path/to/scaffold/directory>`


## User Editable Variables
- `binDir` - Location of software executables, eg: /home/kwade/bin
- `projectDir` - Loctation to send MHConstructor output, eg: 
/home/kwade/MHC_wgs
- `targetCapture` - Path to target-capture MHC fastq, eg: 
/home/TargetCaptureFastq, if none, use 0
- `WGSdata` - Run on local, WGS sequence data in .cram format, yes (1), no (0)
- `oneKGPdata` - To run on 1000genomeProject FTP .cram files, provide 
path to file containing .cram file ftp links, eg: 
$/home/MHC_wgs/1000genomesFtps.txt, if not desired, set to 0
- `HLAgenotypes` - Path to .csv/tab delim .txt containing HLADRB genotypes
- `C4genotypes` - Path to .csv/tab delim .txt containing C4 genotypes
- `R1ext` - File extension to locate target-capture fastq files, eg: R1_001.fastq.gz,
if none, use 0
- `R2ext` - File extension to locate target-capture fastq files, eg: R2_001.fastq.gz,
if none, use 0

- `chromPos` - Coordinates of desired MHC region on hg38, eg; 
"chr6:28509120-33481577"
- `nThreads` - Number of parallel multi threads to use
- `insertSize` - Average mapped read insert size (bp), eg: 400
- `expCov` - Expected average coverage depth across region, eg: 30
- `kmer` - Kmer size for de novo assembly, eg: 51
- `kmerScaff` - Kmer size for de novo scaffolding, eg: 61
- `readCount` - Number of starting reads (in millions) to subsample fastqs, eg: 
5, or, if no subsampling required, set to ‘all’
- `assignHaps` - Run script to assign BMHs, yes (1), no (0)
