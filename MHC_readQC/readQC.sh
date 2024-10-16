#!/bin/bash
## Read QC
## Poly-G removal, read subsampling (if desired), and QC filtering
## Edited: 9/11/24
## Authors: wadekj and rsuseno

### Input arguments ###
# $1 - full R1.fastq file location
# $2 - R1 ID name
# $3 - full R2.fastq file location
# $4 - R2 ID name
# $5 - workPathDir
# $6 - n threads
# $7 - readCount (in millions)
# $8 - sample ID
# $9 - fastqDir

progFastQC=fastqc ##conda
progTrimmomatic=trimmomatic ##conda
readLog=${5}/log_${8}_readQC.txt
million=1000000
declare -i readNum=$7
readInt=`expr ${readNum} \* ${million}`
readIntHalf=`expr ${readInt} / 2`


## Remove stretches of poly-g artifacts ####
fastp -i ${1} -o ${2}_pg.fastq.gz -I ${3} -O ${4}_pg.fastq.gz --trim_poly_g --poly_g_min_len 20 

## Check for number of reads and filter if needed ####
declare -i readNum=$(zless ${2}_pg.fastq.gz | grep '+' | wc -l)
# readNum=$(zless "${2}_pg.fastq.gz" | grep '+' | wc -l)
if [ ${readNum} -ge ${readInt} ]
then	
	seqtk sample -s 347 ${2}_pg.fastq.gz ${readIntHalf} > ${2}_pg_${7}.fastq
	seqtk sample -s 347 ${4}_pg.fastq.gz ${readIntHalf} > ${4}_pg_${7}.fastq
fi

if [ ${readNum} -lt ${readInt} ]
then
	mv ${2}_pg.fastq.gz ${2}_pg_${7}.fastq.gz
	mv ${4}_pg.fastq.gz ${4}_pg_${7}.fastq.gz
fi
gzip ${2}_pg_${7}.fastq
gzip ${4}_pg_${7}.fastq


names=(${8})
reads1=(${2}_pg_${7}.fastq.gz)
reads2=(${4}_pg_${7}.fastq.gz)

cd ${9}


## Adapted from:
#########################################################
#  Reference-guided de novo assembly - SOAP
# ====================================================
# by Heidi Lischer, 2015/2016
# https://bitbucket.org/HeidiLischer/refguideddenovoassembly_pipelines/src/master/
#########################################################

# 1. Step: quality/adapter trimming and quality check:
#######################################################
  # quality check ----------
  echo "quality check of raw reads..."
  for i in ${!names[*]};  #for all indexes in the array
  do
    ${progFastQC} -t ${6} -o ${5} ${reads1[i]} ${reads2[i]}
  done
  
  # quality/adapter trimming --------
  # - remove Illumina adapters provided in the primer files
  # - remove leading and trailing low quality basses (<3) or N
  # - 4 base sliding window -> remove when average quality is < 15
  # - remove reads which are shorter than 40 bp
  echo "quality/adapter trimming..." 
  cd $5
  read1TrimPair=()
  read1TrimUnPair=()
  read2TrimPair=()
  read2TrimUnPair=()
  for i in ${!names[*]}  #for all indexes in the array
  do
    read1TrimPair[i]=${5}/${names[i]}_R1_trimPair.fastq
    read1TrimUnPair[i]=${5}/${names[i]}_R1_trimUnPair.fastq
    read2TrimPair[i]=${5}/${names[i]}_R2_trimPair.fastq
    read2TrimUnPair[i]=${5}/${names[i]}_R2_trimUnPair.fastq

    ${progTrimmomatic} PE -threads ${6} ${reads1[i]} ${reads2[i]} ${read1TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimPair[i]} ${read2TrimUnPair[i]} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 

  done


  # quality check ----------
  echo "quality check of trimmed reads..."
  for i in ${!names[*]}  #for all indexes in the array
  do
    ${progFastQC} -t ${6} -o ${5} ${read1TrimPair[i]} ${read2TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimUnPair[i]}
  done


