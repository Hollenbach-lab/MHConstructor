# MHConstructor

## Project Description

MHConstructor is an high-throughput, haplotype-informed, MHC-specific, short read de novo assembly method intended to further understand the role of MHC in human traits and disease studies. The full article of [MHConstructor can be found here.](https://stock.adobe.com/search?k=cat)

## Installation Guide
To obtain the project directory, the repository can be cloned using the following command:

```
git clone https://github.com/wadekj/MHConstructor.git
cd ./MHConstructor
```

To ensure a reliable run, we have containerized our image in [Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/installation.html). We have tested it on Singularity version 3.11.4. The Singularity recipe file, `MHConstructor/container/mhconstructor.def`, can be built using one of the following commands:
```
cd ./container
singularity build --fakeroot mhconstructor.sif mhconstructor.def
cd ../
```


This command can take approximately an hour to create a .sif file that will need to be included when running the pipeline through Singularity.


## Usage Guide
### Updating user-edtiable variables
To specify various variables for your run, you can edit the `MHConstructor/control.txt` file to your specific need. This is where you specify information such as the location of your data, the location of the output, and many more as described on the [User Editable Variables](#user-editable-variables) section below.

### Run end-to-end pipeline
After you specified the inputs, you can run the following line from within the MHConstructor directory to execute the whole pipeline:

`singularity exec --bind <fastq_location> container/mhconstructor.sif /bin/bash MHCgenerate_v5.sh testID.txt`

### Run specific submodules of interest
In case you are only interested in a specific submodule for your research or analysis purpose, each of the process can be run separately using the following lines:
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
