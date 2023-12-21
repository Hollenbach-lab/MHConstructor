#!/bin/bash
## Wrapper script to contain the full workflow of MHC assembly generation, var calls
## C4 genotyping, HLA genotyping, reference haplotype assigning, 
## pseudo-phasing, non-HLA gene calls, and association
## 2/05/23
## Revised v3: 08/09/23
## wadekj

## Usage sh generateMHC.sh [1] idsFile.txt

###########################################################

########### User editable variables ###########
### Converted to control.txt separate file ###
source ./control.txt


##### DO NOT CHANGE BELOW THIS LINE ########################################
repoDir=${binDir}/MHCconstructor
progSamtools=${repoDir}/tools/samtools-1.9/samtools
progBowtie2=${repoDir}/tools/bowtie2-2.4.2-sra-linux-x86_64/bowtie2
progPicard=picard
assemblyDir=${projectDir}/MHC_assembly
consDir=${projectDir}/consensusHap
samDir=${projectDir}/sam
RefMHChaplotypes=${repoDir}/MHC_hapAssign/mhcRefHaps.txt
mkdir ${projectDir}
mkdir ${assemblyDir}
mkdir ${consDir}
mkdir ${samDir}
conda activate amosPy27
conda env list


if [ $assignHaps == 0 ]
then
	$Hap1== 'hg38'
	$Hap2== 'hg38'
fi
if [ $assignHaps == 1 ]
then
	mkdir ${binDir}/MHC_generate/MHC_hapAssign/bestHaps
	cd MHC_hapAssign
	python assignMHChaps.py ${RefMHChaplotypes} ${HLAgenotypes} ${C4genotypes}
	cd ..
fi

#conda activate amosPy27
#conda env list
ids="$(cat $1)"
for i in $ids; do
	echo $i
        sampleAssembly=${assemblyDir}/${i}_athena_v3
        mkdir $sampleAssembly


######### Pre-process: Step 0.A. Check which type of data the user has, extract WGS reads if needed ######
	## Option A: In-house WGS data that has been aligned to a reference genome (.cram)
	if [ "${WGSdata}" !=0 ]
	then
		cd ./MHC_wgsPreProcess
		sh extractCram.sh ${i} ${chromPos} ${progSamtools}
                R1=${i}_mhc_all_R1.fastq.gz
                R2=${i}_mhc_all_R2.fastq.gz
                R1A="$(echo $R1 | cut -d '.' -f 1,1)"
                R2A="$(echo $R2 | cut -d '.' -f 1,1)"
		fastqDir=${sampleAssembly}
		cd ..
	fi

	## Option B: 1000GenomesProject data via ftp server (list of ftps)
	if [ "${OneKGP}" !=0 ]
	then
		cd ./MHC_wgsPreProcess
		sh extractFrom1KGP.sh ${i} ${oneKGPdata} ${chromPos} ${progSamtools} ${sampleAssembly} 
		R1=${i}_mhc_all_R1.fastq.gz
		R2=${i}_mhc_all_R2.fastq.gz
                R1A="$(echo $R1 | cut -d '.' -f 1,1)"
                R2A="$(echo $R2 | cut -d '.' -f 1,1)"
		fastqDir=${sampleAssembly}
		cd ..
	fi
	
	## Option C: in-house target capture data as .fastq.gz ###
	if [ "${targetCapture}" != 0 ]
	then
                R1="$(find ${fastqDir} -name ${i}*${R1ext})"
                R2="$(find ${fastqDir} -name ${i}*${R2ext})"
		R1A="$(echo $R1 | cut -d '.' -f 1,1)"
		R2A="$(echo $R2 | cut -d '.' -f 1,1)"
		fastqDir=${targetCapture}

############ Preprocess Step 0.B: Assign closest 2 reference MHC haplotypes ###########
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

############## 1. Read quality filter ###########
         cd ./MHC_readQC
         readLog=$sampleAssembly/log_${i}_readQC.txt
         sh readQC.sh ${R1} ${R1A} ${R2} ${R2A} ${sampleAssembly} ${nThreads} ${readCount} ${i} ${fastqDir} >> $readLog
         cd ..

############## 2. Generate ref guided, de novo MHC assembly ##########
	## One assembly for each MHC haplotype. 
	## If heterozygous, two assemblies
	## If homozygous, one assembly
 
	## Generate assembly for Hap1
	cd ./MHC_assembly/
	assemblyLog1=${sampleAssembly}/log_${i}_${Hap1}_assembly.txt
	sh mhc_assembly.sh ${i} ${R1A}_pg_${readCount}.fastq.gz ${R2A}_pg_${readCount}.fastq.gz ${readCount} ${fastqDir} ${Hap1} ${sampleAssembly} ${binDir} ${progSamtools} ${progPicard} ${insertSize} ${kmerAssembly} ${expCov} ${Hap2} ${repoDir} ${progBowtie2}>> $assemblyLog1 
	
	## Check if individual is MHC homozygous or heterozygous
	## If heterozygous, repeat assembly using Hap2 as guide 
	if [ ${Hap1} != ${Hap2} ]
	then
		assemblyLog2=log_${i}_${Hap2}_assembly.txt
		sh mhc_assembly.sh ${i} ${R1A}_pg_${readCount}.fastq.gz ${R2A}_pg_${readCount}.fastq.gz ${readCount} ${fastqDir} ${Hap2} ${sampleAssembly} ${binDir} ${progSamtools} ${progPicard} ${insertSize} ${kmerAssembly} ${expCov} ${Hap1} ${repoDir} ${progBowtie2} >> $assemblyLog2
	
	fi
	cd ..

	####### 3. Error-correct and scaffold assembly ##########
	cd ./MHC_scaffold/
	sh mhc_scaffold.sh $i ${R1A} ${R2A} ${fastqDir} ${sampleAssembly} ${insertSize} ${repoDir} ${Hap1} ${Hap2} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount} 
        
	if [ ${Hap1} != ${Hap2} ]
	then

		sh mhc_scaffold.sh $i ${R1A} ${R2A} ${fastqDir} ${sampleAssembly} ${insertSize} ${repoDir} ${Hap2} ${Hap1} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount}
	fi
	cd ..

	source activate py35
	conda env list

	###### 4. Scaffold and order assembly against refHap(s) #######
	cd ./MHC_order/
        sh orderAssembly.sh ${sampleAssembly}/${i}_${readCount}_${Hap1}/merged_corr/scaffold_gapClosed/${i}_63_final.fasta ${Hap1} ${sampleAssembly}/${i}_${readCount}_${Hap1}/merged_corr/scaffold_gapClosed
        mv ${sampleAssembly}/${i}_${readCount}_${Hap1}/merged_corr/scaffold_gapClosed/RT_${Hap1}/ragtag.scaffold.fasta ${projectDir}/consensusHap/${i}_vs_${Hap1}.fasta

        if [ ${Hap1} != ${Hap2} ]
        then
		sh orderAssembly.sh ${sampleAssembly}/${i}_${readCount}_${Hap2}/merged_corr/scaffold_gapClosed/${i}_63_final.fasta ${Hap2} ${sampleAssembly}/${i}_${readCount}_${Hap2}/merged_corr/scaffold_gapClosed
         	mv ${sampleAssembly}/${i}_${readCount}_${Hap2}/merged_corr/scaffold_gapClosed/RT_${Hap2}/ragtag.scaffold.fasta ${projectDir}/consensusHap/${i}_vs_${Hap2}.fasta
	fi

	cd ..

