#!/bin/bash
## Wrapper script to contain the full workflow of MHC assembly generation, var calls
## C4 genotyping, HLA genotyping, reference haplotype assigning, 
## pseudo-phasing, non-HLA gene calls, and association
## 2/05/23
## Revised v3: 08/09/23
## wadekj

## Usage sh generateMHC.sh [1] idsFile.txt

###########################################################
o
########### User editable variables ###########
## Convert to separate control file, read-in variables? ####
binDir=/home/kwade/bin
projectDir=/home/kwade/MHC_wgs
WGSdir='0'
#oneKGPdata=${projectDir}/ftpLinks_101223.txt
oneKGPdata=${binDir}/MHC_generate_wgs/ftpNA20129.txt
chromPos="chr6:28509120-33481577"
#chromPos="chr6:29500000-33500000"
nThreads=8
insertSize=400
expCov=30
kmer=51
primerFile=test.fasta
readCount=all
assignHaps=1
##**HLAgenotypes=/home/kwade/3_MSgenosOmixon/completedHLA_C4paired_23feb23.csv
##**C4genotypes=/home/kwade/3_MSgenosC4/completedC4_hlapaired_23feb23.csv
progSamtools=/usr/bin/samtools-1.9/samtools
#progPicard=/home/picard/picard/build/libs/picard.jar
progPicard=picard
progBowtie2=/home/kwade/bin/bowtie2-2.4.2-sra-linux-x86_64/bowtie2

##### DO NOT CHANGE BELOW THIS LINE ########################################
#fastqDir=${projectDir}/fastqs
assemblyDir=${projectDir}/4_assemblyAthenav3
denovoVCDir=${projectDir}/5_denovoVarCallv3
consDir=${denovoVCDir}/consensusHap
pafDir=${denovoVCDir}/paf
samDir=${denovoVCDir}/sam
vcfDir=${denovoVCDir}/vcf
hetCalls=${denovoVCDir}/hetCalls
indel=${denovoVCDir}/indel
RefMHChaplotypes=${binDir}/MHC_generate/MHC_hapAssign/mhcRefHaps.txt
#mkdir ${fastqDir}
mkdir ${assemblyDir}
mkdir ${denovoVCDir}
mkdir ${consDir}
mkdir ${pafDir}
mkdir ${samDir}
mkdir ${vcfDir}
mkdir ${hetCalls}
mkdir ${indel}
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
        ######## 1. Assign closest 2 reference MHC haplotypes ###########
        ## ** Requires user to have HLA and C4 genotype files .txt **
        ## A) mHC_hapAssign(C4 genotypes, HLA genotypes) --> .txt
        if [ ${assignHaps} -eq 1 ]
        then
                cd ./MHC_hapAssign/
                alleles="$(cat ${binDir}/MHC_generate_wgs/MHC_hapAssign/bestHaps/${i}_bestMHChaps.txt)";
                Hap1="$(echo $alleles| cut -d ' ' -f 2,2)"
                Hap2="$(echo $alleles| cut -d ' ' -f 5,5)"
                echo $Hap1
                echo $Hap2
                cd ..
        fi

######### Pre-process: Check which type of data the uder has ######
	
	## A: In-house WGS data that has been aligned to a reference genome (.cram)
	if [ "$WGSdata" != '0' ]
	then
		cd ./MHC_wgsPreProcess
		sh extractCram.sh ${i} ${chromPos} ${progSamtools}
		cd ..
	fi

	## B: 1000GenomesProject data via ftp server (list of ftps)
	if [ "${OneKGP}" != '0' ]
	then
		
		cd ./MHC_wgsPreProcess
		sh extractFrom1KGP.sh ${i} ${oneKGPdata} ${chromPos} ${progSamtools} ${sampleAssembly} 
		cd ..
	fi
	
	######## 0. Read quality filter ###########
	 cd ./MHC_readQC
	 readLog=$sampleAssembly/log_${i}_readQC.txt
	 /bin/bash readQC.sh ${i} ${sampleAssembly}/${i}_mhc_all_R1.fastq ${sampleAssembly}/${i}_mhc_all_R2.fastq ${sampleAssembly} ${nThreads} >> $readLog
	 cd ..
#########################################################################

	######## 2. Generate ref guided, de novo MHC assembly ##########
	## One assembly for each MHC haplotype. 
	## If heterozygous, two assemblies
	## If homozygous, one assembly
 
	## Generate assembly for Hap1
	cd ./MHC_assembly23/
	echo 'MHC_assembly23'
	echo $i
	assemblyLog1=${sampleAssembly}/log_${i}_${Hap1}_assembly.txt
	sh run_athena_v3.sh ${i} ${sampleAssembly}/${i}_mhc_all_R1 ${sampleAssembly}/${i}_mhc_all_R2 ${sampleAssembly} ${Hap1} ${sampleAssembly} ${binDir} ${progSamtools} ${progPicard} ${readCount} ${insertSize} ${kmer} ${expCov} ${Hap2} ${binDir} ${progBowtie2}>> $assemblyLog1 
	
	## Check if individual is MHC homozygous or heterozygous
	## If heterozygous, repeat assembly using Hap2 as guide 
	if [ ${Hap1} != ${Hap2} ]
	then
		assemblyLog2=log_${i}_${Hap2}_assembly.txt
		sh run_athena_v3.sh ${i}  ${sampleAssembly}/${i}_mhc_all_R1 ${sampleAssembly}/${i}_mhc_all_R2 ${sampleAssembly} ${Hap2} ${sampleAssembly} ${binDir} ${progSamtools} ${progPicard} ${readCount} ${insertSize} ${kmer} ${expCov} ${Hap1} ${binDir} ${progBowtie2} >> $assemblyLog2
	
	fi
	cd ..

	####### 3. Error-correct and scaffold assembly ##########
	cd ./MHC_scaffold/
	sh refGuidedDeNovoAssembly_velvet_athena_scaffold_v3_target_haploid.sh $i ${sampleAssembly}/${i}_mhc_all_R1 ${sampleAssembly}/${i}_mhc_all_R2 ${sampleAssembly} ${sampleAssembly} ${insertSize} ${binDir} ${Hap1} ${Hap2} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount} 
        
	if [ ${Hap1} != ${Hap2} ]
	then

		sh refGuidedDeNovoAssembly_velvet_athena_scaffold_v3_target_haploid.sh $i ${sampleAssembly}/${i}_mhc_all_R1 ${sampleAssembly}/${i}_mhc_all_R2 ${sampleAssembly} ${sampleAssembly} ${insertSize} ${binDir} ${Hap2} ${Hap1} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount}
	fi
	cd ..

	#conda activate py35

	###### 4. Scaffold and order assembly against refHap(s) #######
	'''cd ./MHC_order/
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
	cd ..'''


	######### 6. Merge scaffolded assembly vs hg38 var calls into 1 .vcf if heterozygous #############
	##cd ./MHC_varHg38
	##sh mergeHg38Var.sh 
	##if [ ${Hap1} != ${Hap2} ]
        ##then
        	##sh mergeHg38Var.sh
	##fi
	##cd ..
	


done
