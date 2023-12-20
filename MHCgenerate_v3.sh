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
## Convert to separate control file, read-in variables? ####
binDir=/home/kwade/bin
projectDir=/home/kwade/MHC_assemblies23_v3
fastqDir=/home/kwade/INDIGO_Fastq
#fastqDir=/home/kwade/normanSRA
#R1ext=R1_001.fastq.gz
#R2ext=R2_001.fastq.gz
R1ext=R1_001.fastq.collapsed.fq.gz
R2ext=R2_001.fastq.collapsed.fq.gz
nThreads=8
insertSize=400
expCov=60
kmer=51
primerFile=test.fasta
readCount=5
assignHaps=1
##**HLAgenotypes=/home/kwade/3_MSgenosOmixon/completedHLA_C4paired_23feb23.csv
##**C4genotypes=/home/kwade/3_MSgenosC4/completedC4_hlapaired_23feb23.csv
progSamtools=/usr/bin/samtools-1.9/samtools
#progPicard=/home/picard/picard/build/libs/picard.jar
progPicard=picard
progBowtie2=/home/kwade/bin/bowtie2-2.4.2-sra-linux-x86_64/bowtie2

##### DO NOT CHANGE BELOW THIS LINE ########################################
assemblyDir=${projectDir}/4_assemblyAthenav3
denovoVCDir=${projectDir}/5_denovoVarCallv3
consDir=${denovoVCDir}/consensusHap
pafDir=${denovoVCDir}/paf
samDir=${denovoVCDir}/sam
vcfDir=${denovoVCDir}/vcf
hetCalls=${denovoVCDir}/hetCalls
RefMHChaplotypes=${binDir}/MHC_generate/MHC_hapAssign/mhcRefHaps.txt
mkdir ${assemblyDir}
mkdir ${denovoVCDir}
mkdir ${consDir}
mkdir ${pafDir}
mkdir ${samDir}
mkdir ${vcfDir}
mkdir ${hetCalls}

## If assignHaps== ###
# if [ $assignHaps == 0 ]
# $Hap1== 'hg38'
# $Hap2== 'hg38'
# fi
## if [ $assignHaps == 1 ]
## mkdir ${binDir}/MHC_generate/MHC_hapAssign/bestHaps
## cd MHC_hapAssign
## echo $$ > parent_pid.txt
#python assignMHChaps.py ${RefMHChaplotypes} ${HLAgenotypes} ${C4genotypes}
## cd ..
# fi

ids="$(cat $1)"
for i in $ids; do
	echo $i
	sampleAssembly=${assemblyDir}/${i}_athena_v3
    	mkdir $sampleAssembly
	
	R1="$(find ${fastqDir} -name ${i}*${R1ext})"
	R2="$(find ${fastqDir} -name ${i}*${R2ext})"
	R1A="$(echo $R1 | cut -d '.' -f 1,1)"
	R2A="$(echo $R2 | cut -d '.' -f 1,1)"
	echo $sampleAssembly
	echo $assemblyDir
	######## 0. Read quality filter ###########
	 cd ./MHC_readQC
	 readLog=$sampleAssembly/log_${i}_readQC.txt
	 /bin/bash readQC.sh ${R1} ${R1A} ${R2} ${R2A} ${sampleAssembly} ${nThreads} ${primerFile} ${readCount} ${i} ${fastqDir} >> $readLog
	 cd ..
	 start
#########################################################################
	######## 1. Assign closest 2 reference MHC haplotypes ###########
	## ** Requires user to have HLA and C4 genotype files .txt **
        ## A) mHC_hapAssign(C4 genotypes, HLA genotypes) --> .txt
	if [ ${assignHaps} -eq 1 ]
	then
		cd ./MHC_hapAssign/
        	alleles="$(cat ${binDir}/MHC_generate_v3/MHC_hapAssign/bestHaps/${i}_bestMHChaps.txt)";
        	Hap1="$(echo $alleles| cut -d ' ' -f 2,2)"
        	Hap2="$(echo $alleles| cut -d ' ' -f 8,8)"
		echo $Hap1
        	echo $Hap2
        	cd ..
	fi 
	start=`date +%s`
	######## 2. Generate ref guided, de novo MHC assembly ##########
	## One assembly for each MHC haplotype. 
	## If heterozygous, two assemblies
	## If homozygous, one assembly
 
	## Generate assembly for Hap1
	cd ./MHC_assembly23/
	echo 'MHC_assembly23'
	echo $i
	assemblyLog1=${sampleAssembly}/log_${i}_${Hap1}_assembly.txt
	sh run_athena_v3.sh ${i} ${R1A} ${R2A} ${fastqDir} ${Hap1} ${sampleAssembly} ${binDir} ${progSamtools} ${progPicard} ${readCount} ${insertSize} ${kmer} ${expCov} ${Hap2} ${binDir} ${progBowtie2}>> $assemblyLog1 
	
	## Check if individual is MHC homozygous or heterozygous
	## If heterozygous, repeat assembly using Hap2 as guide 
	if [ ${Hap1} != ${Hap2} ]
	then
		assemblyLog2=log_${i}_${Hap2}_assembly.txt
		sh run_athena_v3.sh ${i} ${R1A} ${R2A} ${fastqDir} ${Hap2} ${sampleAssembly} ${binDir} ${progSamtools} ${progPicard} ${readCount} ${insertSize} ${kmer} ${expCov} ${Hap1} ${binDir} ${progBowtie2} >> $assemblyLog2
	
	fi
	cd ..
	end=`date +%s`
	runtime=$((end-start))
	echo 'Assembly runtime:'
	echo $runtime

        start=`date +%s`

	####### 3. Error-correct and scaffold assembly ##########
	cd ./MHC_scaffold/
	sh refGuidedDeNovoAssembly_velvet_athena_scaffold_v3_target_haploid.sh $i ${R1A} ${R2A} ${fastqDir} ${sampleAssembly} ${insertSize} ${binDir} ${Hap1} ${Hap2} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount} 
        
	if [ ${Hap1} != ${Hap2} ]
	then

		sh refGuidedDeNovoAssembly_velvet_athena_scaffold_v3_target_haploid.sh $i ${R1A} ${R2A} ${fastqDir} ${sampleAssembly} ${insertSize} ${binDir} ${Hap2} ${Hap1} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount}
	fi
	cd ..
	end=`date +%s`
        runtime=$((end-start))
        echo 'Scaffold runtime:'
        echo $runtime
	
	#conda activate py35
	echo $Hap1
	echo $Hap2	
	###### 4. Scaffold and order assembly against refHap(s) #######
	'''cd ./MHC_order/
        ## add source instead of sh to run py35 env
	sh orderAssembly.sh ${sampleAssembly}/${i}_${readCount}_${Hap1}/merged_corr/scaffold_gapClosed/${i}_63_final.fasta ${Hap1} ${sampleAssembly}/${i}_${readCount}_${Hap1}/merged_corr/scaffold_gapClosed
        mv ${sampleAssembly}/${i}_${readCount}_${Hap1}/merged_corr/scaffold_gapClosed/RT_${Hap1}/ragtag.scaffold.fasta ${denovoVCDir}/consensusHap/${i}_vs_${Hap1}.fasta

        if [ ${Hap1} != ${Hap2} ]
        then
		sh orderAssembly.sh ${sampleAssembly}/${i}_${readCount}_${Hap2}/merged_corr/scaffold_gapClosed/${i}_63_final.fasta ${Hap2} ${sampleAssembly}/${i}_${readCount}_${Hap2}/merged_corr/scaffold_gapClosed
         	mv ${sampleAssembly}/${i}_${readCount}_${Hap2}/merged_corr/scaffold_gapClosed/RT_${Hap2}/ragtag.scaffold.fasta ${denovoVCDir}/consensusHap/${i}_vs_${Hap2}.fasta
	fi

	cd ..

	######### 5. Generate denovo varCalls ##########################
	## A) denovovarcall.sh(MHC hap1, MHC hap2, assembly.fasta) --> dnHap1.vcf, dnHap2.vcf 
	cd ./MHC_denovoVarCall
	sh run_denovo_varCall_samtools.sh ${i} ${Hap1} ${denovoVCDir}/consensusHap/${i}_vs_${Hap1}.fasta ${denovoVCDir} ${insertSize} ${sampleAssembly} 2> ${denovoVCDir}/indel/${i}_vs_${Hap1}_dn_indels.txt	
        
	if [ ${Hap1} != ${Hap2} ]
        then
                sh run_denovo_varCall_samtools.sh ${i} ${Hap2} ${denovoVCDir}/consensusHap/${i}_vs_${Hap2}.fasta ${denovoVCDir} ${insertSize} ${sampleAssembly} 2> ${denovoVCDir}/indel/${i}_vs_${Hap2}_dn_indels.txt
	fi
	cd ..


	######### 6. Merge scaffolded assembly vs hg38 var calls into 1 .vcf if heterozygous #############
	##cd ./MHC_varHg38
	##sh mergeHg38Var.sh 
	##if [ ${Hap1} != ${Hap2} ]
        ##then
        	##sh mergeHg38Var.sh
	##fi
	##cd ..'''
	


done
