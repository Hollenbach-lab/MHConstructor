#! /bin/bash
## Script to access the ftp location of a .cram file and extract MHC
## specifically designed to access 1000genomesProject data
## 09/12/23
## wadekj

# $1 - sample ID
# $2 - ftp location list file
# $3 - chromosome position, default: chr6:28525013-33457522
# $4 - samtools location
# $5 - sampleAssembly

progSamtools=${4}
altChroms="$(cat GRCH38chroms_altMHC.txt)"
#ftp=$(python getFTP.py ${1} ${2} 2>&1)

echo 'Downloading reference fasta'
wget -nc https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
echo 'Finished downloading reference fasta'

# Extract the ftp link associated with a sample

ftp=$(less ${2} | grep ${1} )  
echo $ftp
mkdir -p ${5}/extractedReads
${progSamtools} view -f 1 -h -T ./GRCh38_full_analysis_set_plus_decoy_hla.fa -b ${ftp} ${3}| ${progSamtools} sort - -n -o ${5}/extractedReads/${1}_mhc_paired.sorted.bam

${progSamtools} fastq -n ${5}/extractedReads/${1}_mhc_paired.sorted.bam -1 ${5}/extractedReads/${1}_mhc_paired_R1.fastq -2 ${5}/extractedReads/${1}_mhc_paired_R2.fastq -0 /dev/null -s /dev/null

# Get alt locus mapped reads
for i in ${altChroms};
do
	altLoc="$(echo $i | cut -d '>' -f 2,2)"
	echo ${altLoc}
	${progSamtools} view -f 1 -h -T ./GRCh38_full_analysis_set_plus_decoy_hla.fa -b ${ftp} ${altLoc} | ${progSamtools} sort - -n -o ${5}/extractedReads/${1}_mhc_${altLoc}_paired.sorted.bam
	${progSamtools} fastq -n ${5}/extractedReads/${1}_mhc_${altLoc}_paired.sorted.bam -1 ${5}/extractedReads/${1}_mhc_${altLoc}_paired_R1.fastq -2 ${5}/extractedReads/${1}_mhc_${altLoc}_paired_R2.fastq -0 /dev/null -s /dev/null
done;

## Get unmapped reads.bam
${progSamtools} view -h -T ./GRCh38_full_analysis_set_plus_decoy_hla.fa -b ${ftp} '*' | ${progSamtools} sort - -n -o ${5}/extractedReads/${1}_unmap.sorted.bam

mv *.crai ${5}

## Unmapped reads
${progSamtools} fastq -n ${5}/extractedReads/${1}_unmap.sorted.bam -1 ${5}/extractedReads/${1}_mhc_unmap_R1.fastq -2 ${5}/extractedReads/${1}_mhc_unmap_R2.fastq -0 /dev/null -s /dev/null


### Merge all extracted reads together
for r1 in ${5}/extractedReads/*R1.fastq;
do 
	cat $r1 >> ${5}/${1}_mhc_all_R1.fastq
done;

for r2 in ${5}/extractedReads/*R2.fastq;
do 
	cat $r2 >>${5}/${1}_mhc_all_R2.fastq
done;

mv ${5}/extractedReads/${1}_mhc_unmap_R1.fastq ${5}
mv ${5}/extractedReads/${1}_mhc_unmap_R2.fastq ${5}

gzip ${5}/${1}_mhc_unmap_R1.fastq
gzip ${5}/${1}_mhc_all_R1.fastq

gzip ${5}/${1}_mhc_unmap_R2.fastq
gzip ${5}/${1}_mhc_all_R2.fastq
