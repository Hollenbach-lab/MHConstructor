# MHConstructor

## Project Description

MHConstructor is an high-throughput, haplotype-informed, MHC-specific, short read de novo assembly method intended to further understand the role of MHC in human traits and disease studies. The full article of [MHConstructor can be found here.](https://www.biorxiv.org/content/10.1101/2024.05.20.595060v1)

## Installation Guide
To obtain the project directory, the repository can be cloned using the following command:

```
git clone https://github.com/Hollenbach-lab/MHConstructor.git
cd ./MHConstructor
```

To ensure a reliable run, we have containerized our image in [Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/installation.html). We have tested it on Singularity version 3.11.4. To obtain the image, you can use one of the following commands:  
1. Build with sudo (50 minutes)
```
cd ./container
sudo singularity build mhconstructor.sif mhconstructor.def
cd ../
```

2. Build with `fakeroot`, if sudo is not possible (50 minutes)
```
cd ./container
singularity build --fakeroot mhconstructor.sif mhconstructor.def
cd ../
```

3. Pull the image directly from Sylabs (25 minutes, depends on your internet connection)
```
cd ./container
singularity pull mhconstructor.sif library://rsuseno/rsuseno/mhconstructor:latest
cd ../
```
An `mhconstructor.sif` file will be created and ready to use!


## Usage Guide
### Updating user-edtiable variables
To specify various variables for your run, you can edit the `MHConstructor/control.txt` file to your specific need. This is where you specify information such as the location of your data, the location of the output, and many more as described on the [User Editable Variables](#user-editable-variables) section below.

### Run end-to-end pipeline
After you specified the inputs, you can run the following line from within the MHConstructor directory to execute the whole pipeline:

```
singularity exec --bind <fastq_location> container/mhconstructor.sif /bin/bash MHCgenerate_v5.sh testID.txt
```

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


## Output
The primary use of this pipeline is to process raw sequencing data and produce the consensus haplotype FASTA. This file can be found under the `consensusHap` directory in the output folder, and will follow a format of `<sample>_vs_<haplotype>.fasta` as its filename. If a sample is heterozygous, it can be mapped to two haplotypes, in which case there would exist two FASTA files for one sample.

Other information such as downsampled reads or assembly files are stored in the `MHConstructor_assemblies` folder, should you be interested in looking into it.
```
├── consensusHap
│   └── EPIC0345_vs_OK649231.fasta
│   └── EPIC0345_vs_OK649232.fasta
├── MHConstructor_assemblies
│   └── EPIC0345_MHConstructor
└── sam
```
