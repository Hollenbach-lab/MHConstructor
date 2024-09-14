# MHConstructor

## Project Description

MHConstructor is an high-throughput, haplotype-informed, MHC-specific, short read de novo assembly method intended to further understand the role of MHC in human traits and disease studies. The full article for [MHConstructor can be found here.](https://www.biorxiv.org/content/10.1101/2024.05.20.595060v1)

## Installation Guide
To obtain the project directory, the repository can be cloned using the following command:

```
git clone https://github.com/Hollenbach-lab/MHConstructor.git
cd ./MHConstructor
```

To ensure a reliable run, we have containerized our image in Singularity (tested on version 3.11.4). You can check whether Singularity has been installed in your system by running `singularity --version`. Please install Singularity by following their guide [here](https://docs.sylabs.io/guides/3.11/admin-guide/installation.html). 

To obtain the image, you can use one of the following commands:  
- Option A: Build with `fakeroot` (~50 minutes)
```
cd ./container
singularity build --fakeroot mhconstructor.sif mhconstructor.def
cd ../
```

- Option B: Pull the image directly from Sylabs (~25 minutes, depends on your internet connection)
```
cd ./container
singularity pull mhconstructor.sif library://rsuseno/rsuseno/mhconstructor:latest
cd ../
```
- Option C: Build with sudo (~50 minutes)
```
cd ./container
sudo singularity build mhconstructor.sif mhconstructor.def
cd ../
```
An `mhconstructor.sif` file will be created and ready to use!


## Usage Guide
### Updating user-editable variables
To specify various variables for your run, you can edit the `MHConstructor/control.txt` file to your specific need. This is where you specify information such as the location of your data, the location of the output, and many more as described on the [User Editable Variables](#user-editable-variables) section below.
### HLA-DRB1 and C4A/B genotyping
HLA-DRB1 genotypes are required to determine each individual's guide BMH. If you already have these genotypes, provide the path to a csv or tab delimited file in the "HLAgenotypes=" field in the control file.

If you would like to generate HLA-DRB1 genotypes, you can run the T1K (Song et al., 2023) software from within the MHConstructor container as follows:
```
cd genotypes/
singularity exec ../container/mhconstructor.sif /bin/bash genotypeDRB.sh
cd ..
```
This will write a file called `genotypes/HLAgenotypes.csv`. Afterwards, update the control file "HLAgenotypes=" field with the absolute path to point to this file.

To generate C4A/B, L/S and copy number genotypes, please use [C4Investigator](https://github.com/Hollenbach-lab/C4Investigator) (Marin et al., 2024). The repository has instruction to run C4Investigator.
Afterwards, update the control file "C4genotypes=" field with the absolute path to the resulting `C4Investigator_c4_summary.csv` file.

For examples of acceptable file formats, please see the example_HLAgenotypes.csv and example_C4Investigator_c4_summary.csv files in the ./genotypes directory

### Run end-to-end pipeline
After you specified the inputs, you can run the following line from within the MHConstructor directory to execute the whole pipeline:

```
singularity exec --bind <fastq_location> container/mhconstructor.sif /bin/bash MHCgenerate_v5.sh testID.txt
```
`testID.txt` specifies which sample you are running.

(TODO: script to automatically generate a testID.txt file from a given directory to ensure easier usage of MHConstructor)

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
- `expCovUnmap` - Expected average coverage depth for re-mapping of intially unmapped reads
- `kmerAssembly` - Kmer size for de novo assembly, eg: 51
- `kmerScaffold` - Kmer size for de novo scaffolding, eg: 61
- `readCount` - Number of starting reads (in millions) to subsample fastqs, eg: 
5, or, if no subsampling required, set to ‘all’
- `assignHaps` - Run script to assign BMHs, yes (1), no (0) (If 0, will assign the hg38 reference MHC as default guide to all individuals)


## Output
The primary use of this pipeline is to process raw sequencing data and produce the consensus haplotype FASTA. This file can be found under the `consensusHap` directory in the output folder, and will follow a format of `<sample>_vs_<haplotype>_RT_onlyPlaced.fasta` as its filename. If a sample is heterozygous, it can be mapped to two haplotypes, in which case there would exist two FASTA files for one sample, scaffolded with respect to each BMH guide sequence. These files can be found:
- MHConstructor/consensusHap/SampleID_RT_onlyPlaced.fasta

These files represent the single, continuous sequence comprised of assembly contigs that RagTag (Alonge et al., 2022) was able to place against the BMH guide sequence, ie; consensus sequences. However, for each assembly, there will be contigs that were not able to be placed. These likely contain larger structural variants and putative novel sequences. As there is no straightforward, high throughput way to integrate these into the final consensus sequence, it is recommended that users investigate the additional contigs manually, post-hoc. These unplaced contigs can be found in the file:
- MHConstructor/consensusHap/unplacedContigs/SampleID_BMH_RT_unPlaced.fasta

There is a large amount of intermediate data generated during the assembly process that may be of interest to users for optimization and/or developing additional functionality. For an example results output directory, please visit our repository at the Center for Open Science https://osf.io/46y5k/. Here, you will also find a description of each of these files.
