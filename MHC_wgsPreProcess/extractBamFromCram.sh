#! /bin/bash

## Script to iterate through all *.cram files in a provided directory ($1)
## and extract region associated with the MHC reference coordinates
## for hg38 reference genome
## wadekj 
## 05/15/23

## User input args:
# $1 - name of directory containing cram files. If current directory, use "./"
	## Otherwise, like this: "/home/user/exampleDir/"
# $2 - reference genome. If file on local server, list full directory and file name: "/home/user/exampleDir/refGenome.fasta"
	## Otherwise, if on remote server, input full ftp path

## Dependencies- requires samtools to be installed

progSamtools=/usr/bin/samtools-1.9/samtools

#for f in ${1}*.cram;
#do 
	fName="$(echo $1 | cut -d '.' -f 1,1)"
	echo $1
	echo $fName
	## Get MHC.bam
	#${progSamtools} view ${f} > ${fname}.bam
	#${progSamtools} sort ${f} > ${fName}.sort.cram
	#echo 'done sort'
	${progSamtools} index ${1}
	echo 'done index'
	#${progSamtools} view -T ${2} -b -o ${fName}_mhc.bam ${1} chr6:28525013-33457522
	${progSamtools} view -h -T ${2} -b ${1} 28509120-33481577 | ${progSamtools} sort -T ${fName} -n -o ${fName}_mhc.bam
	## Get C4.bam
        #${progSamtools} view -T ${2} -b -o ${fName}_c4.bam ${f}  chr6:

	## Get unmapped reads.bam
        ${progSamtools} view -h -T ${2} -b -o ${fName}_unmap.bam ${1} '*'

