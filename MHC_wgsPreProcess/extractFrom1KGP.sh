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

#ftp=$(python getFTP.py ${1} ${2} 2>&1)

# Extract the ftp link associated with a sample

ftp=$(less ${2} | grep ${1} )  
echo $ftp
${progSamtools} view -f 1 -h -T ../refMHChaps/GRCh38_full_analysis_set_plus_decoy_hla.fa -b ${ftp} ${3}| ${progSamtools} sort - -n -o ${5}/${1}_mhc_paired.sorted.bam 

## Get unmapped reads.bam
${progSamtools} view -h -T ../refMHChaps/GRCh38_full_analysis_set_plus_decoy_hla.fa -b ${ftp} '*' | ${progSamtools} sort - -n -o ${5}/${1}_unmap.sorted.bam

mv *.crai ${5}

${progSamtools} fastq -n ${5}/${1}_mhc_paired.sorted.bam -1 ${5}/${1}_mhc_paired_R1.fastq -2 ${5}/${1}_mhc_paired_R2.fastq -0 /dev/null -s /dev/null

${progSamtools} fastq -n ${5}/${1}_unmap.sorted.bam -1 ${5}/${1}_mhc_unmap_R1.fastq -2 ${5}/${1}_mhc_unmap_R2.fastq -0 /dev/null -s /dev/null

### Merge the mapped and unmapped reads together
cat ${5}/${1}_mhc_paired_R1.fastq ${5}/${1}_mhc_unmap_R1.fastq >> ${5}/${1}_mhc_all_R1.fastq

cat ${5}/${1}_mhc_paired_R2.fastq ${5}/${1}_mhc_unmap_R2.fastq >> ${5}/${1}_mhc_all_R2.fastq
