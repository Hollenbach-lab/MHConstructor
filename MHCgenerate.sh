#!/bin/bash
## Wrapper script to contain the full workflow of MHC assembly generation, var calls
## Edited: 09/11/24
## Authors: wadekj and susenor

## Usage: sh MHCgenerate.sh [1] idsFile.txt

###########################################################

#### Import user editable variables from control file ###########
source activate
source ./control.txt


##### DO NOT CHANGE THESE REPOSITORIES ########################################
repoDir=${binDir}/MHConstructor
progSamtools=${repoDir}/tools/samtools-1.9/samtools
progBowtie2=${repoDir}/tools/bowtie2-2.4.2-sra-linux-x86_64/bowtie2
progPicard=picard
assemblyDir=${projectDir}/MHConstructor_assemblies
consDir=${projectDir}/consensusHap
unplacedContigs=${projectDir}/consensusHap/unplacedContigs
samDir=${projectDir}/sam
RefMHChaplotypes=${repoDir}/MHC_hapAssign/mhcRefHaps.txt
mkdir ${projectDir}
mkdir ${assemblyDir}
mkdir ${consDir}
mkdir ${unplacedContigs}
mkdir ${samDir}
##########################################################################


## If HLA-DRB1 and C4 genotypes have been created, assign BMH
if [ ${assignHaps} -eq 1 ]
then
	if [ ! -d "${repoDir}/MHC_hapAssign/bestHaps" ]
	then 
		mkdir ${repoDir}/MHC_hapAssign/bestHaps
		cd MHC_hapAssign
		python assignMHChaps.py ${RefMHChaplotypes} ${HLAgenotypes} ${C4genotypes}
		cd ..
	fi
fi

## If no HLA-DRB1 and C4 genotypes, don't assign BMH and use the hg38 MHC reference as default
if [ ${assignHaps} -eq 0 ]
then
	Hap1="hg38"
	Hap2="hg38"
fi


### Iterate through input sample list by ID number and perform MHConstructor functions
ids="$(cat $1)"
for i in $ids; do
	echo $i
	conda activate amosPy27
	conda env list
	sampleAssembly=${assemblyDir}/${i}_MHConstructor
	mkdir $sampleAssembly


### Pre-process: Step 0.A. Check which type of short-read data the user has, extract reads if needed ######
	
## Option A: In-house WGS data that has been aligned to a reference genome (.cram)
	## **Not fully de-bugged yet**
	if [ "${WGSdata}" != 0 ]
	then
		cd ./MHC_wgsPreProcess
		sh extractBamFromCram.sh ${i} ${chromPos} ${progSamtools}
		R1=${i}_mhc_all_R1.fastq.gz
		R2=${i}_mhc_all_R2.fastq.gz
		R1A="$(echo $R1 | cut -d '.' -f 1,1)"
		R2A="$(echo $R2 | cut -d '.' -f 1,1)"
		fastqDir=${sampleAssembly}
		cd ..
	fi
	echo ${oneKGPdata}

## Option B: Extract 1000GenomesProject data via ftp server, accessed with input file containing list of ftp locations
	if [ "$oneKGPdata" != "0" ];
	then
		cd ./MHC_wgsPreProcess
		sh extractFrom1KGP.sh ${i} ${oneKGPdata} ${chromPos} ${progSamtools} ${sampleAssembly} 
		echo 'OneKGP'
		R1=${sampleAssembly}/${i}_mhc_all_R1.fastq.gz
		R2=${sampleAssembly}/${i}_mhc_all_R2.fastq.gz
		R1A="$(echo $R1 | cut -d '.' -f 1,1)"
		R2A="$(echo $R2 | cut -d '.' -f 1,1)"
		fastqDir=${sampleAssembly}
		cd ..
	fi
	
## Option C: In-house target capture short read data as .fastq.gz files ###
	if [ "${targetCapture}" != 0 ]
	then
		fastqDir=${targetCapture}
		R1="$(find ${fastqDir} -name ${i}*${R1ext})"
		R2="$(find ${fastqDir} -name ${i}*${R2ext})"
		R1A="$(echo $R1 | cut -d '.' -f 1,1)"
		R2A="$(echo $R2 | cut -d '.' -f 1,1)"

	fi

#### Preprocess Step 0.B: Asscess bestHap.txt files created above and assign closest 2 reference MHC haplotypes as variables $Hap1 and $Hap2 ###########
	if [ ${assignHaps} -eq 1 ]
	then
		cd ./MHC_hapAssign/
		alleles="$(cat ${repoDir}/MHC_hapAssign/bestHaps/${i}_bestMHChaps.txt)";
		Hap1="$(echo $alleles| cut -d ' ' -f 2,2)"
		Hap2="$(echo $alleles| cut -d ' ' -f 8,8)"
		echo $Hap1
		echo $Hap2
		cd ..
	fi 


#### Step 1. Read quality filter, code adapted from Lischer and Shimizu, 2017 #########
	cd ./MHC_readQC
	readLog=$sampleAssembly/log_${i}_readQC.txt
	sh readQC.sh ${R1} ${R1A} ${R2} ${R2A} ${sampleAssembly} ${nThreads} ${readCount} ${i} ${fastqDir} >> $readLog
	cd ..


#### Step 2. Generate ref guided, de novo MHC assembly. Code adapted from Lischer and Shimizu, 2017 ##########
	## One assembly for each MHC haplotype. 
	## If heterozygous BMHs, two assemblies
	## If homozygous BMH, one assembly
 
	## Generate assembly for Hap1
	cd ./MHC_assembly/
	assemblyLog1=${sampleAssembly}/log_${i}_${Hap1}_assembly.txt
	sh mhc_assembly.sh ${i} ${R1A}_pg_${readCount}.fastq.gz ${R2A}_pg_${readCount}.fastq.gz ${readCount} ${fastqDir} ${Hap1} ${sampleAssembly} ${binDir} ${progSamtools} ${progPicard} ${insertSize} ${kmerAssembly} ${expCov} ${Hap2} ${repoDir} ${progBowtie2} ${expCovUnmap}>> $assemblyLog1 
	
	## Check if individual is MHC homozygous or heterozygous
	## If heterozygous, repeat assembly using Hap2 as guide 
	if [ ${Hap1} != ${Hap2} ]
	then
		assemblyLog2=${sampleAssembly}/log_${i}_${Hap2}_assembly.txt
		sh mhc_assembly.sh ${i} ${R1A}_pg_${readCount}.fastq.gz ${R2A}_pg_${readCount}.fastq.gz ${readCount} ${fastqDir} ${Hap2} ${sampleAssembly} ${binDir} ${progSamtools} ${progPicard} ${insertSize} ${kmerAssembly} ${expCov} ${Hap1} ${repoDir} ${progBowtie2} ${expCovUnmap} >> $assemblyLog2
	
	fi
	cd ..


#### Step 3. Error-correct and scaffold assembly. Code adapted from Lischer and Shimizu, 2017 ##########
	cd ./MHC_scaffold/
	sh mhc_scaffold.sh $i ${R1A} ${R2A} ${fastqDir} ${sampleAssembly} ${insertSize} ${repoDir} ${Hap1} ${Hap2} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount} ${kmerScaffold} 
	
	if [ ${Hap1} != ${Hap2} ]
	then

		sh mhc_scaffold.sh $i ${R1A} ${R2A} ${fastqDir} ${sampleAssembly} ${insertSize} ${repoDir} ${Hap2} ${Hap1} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount} ${kmerScaffold}
	fi
	cd ..

	## Switch conda envs to enable use of RagTag
	conda activate py35
	conda env list


#### Step 4. Scaffold and order assembly against refHap(s). #######
	cd ./MHC_order/
	sh orderAssembly.sh ${sampleAssembly}/${i}_${readCount}_${Hap1}/merged_corr/scaffold_gapClosed/${i}_63_final.fasta ${Hap1} ${sampleAssembly}/${i}_${readCount}_${Hap1}/merged_corr/scaffold_gapClosed ${i} ${repoDir}

	if [ ${Hap1} != ${Hap2} ]
	then
		sh orderAssembly.sh ${sampleAssembly}/${i}_${readCount}_${Hap2}/merged_corr/scaffold_gapClosed/${i}_63_final.fasta ${Hap2} ${sampleAssembly}/${i}_${readCount}_${Hap2}/merged_corr/scaffold_gapClosed ${i} ${repoDir}
	fi

	cd ..

done
