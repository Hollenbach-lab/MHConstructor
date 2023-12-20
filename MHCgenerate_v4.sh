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
# binDir=/home/rsuseno
# projectDir=/home/rsuseno/output_MHC/output_MHC_sing_1205
# fastqDir=/home/kwade/INDIGO_Fastq
# R1ext=R1_001.fastq.gz
# R2ext=R2_001.fastq.gz
# nThreads=8
# insertSize=400
# expCov=60
# kmer=51
# primerFile=test.fasta
# readCount=5
# assignHaps=1
##**HLAgenotypes=/home/kwade/3_MSgenosOmixon/completedHLA_C4paired_23feb23.csv
##**C4genotypes=/home/kwade/3_MSgenosC4/completedC4_hlapaired_23feb23.csv


##### DO NOT CHANGE BELOW THIS LINE ########################################
repoDir=${binDir}/MHC_generate_targetCapture_v4
progSamtools=${repoDir}/tools/samtools-1.9/samtools
progBowtie2=${repoDir}/tools/bowtie2-2.4.2-sra-linux-x86_64/bowtie2
progPicard=picard
assemblyDir=${projectDir}/4_assemblyAthenav3
denovoVCDir=${projectDir}/5_denovoVarCallv3
consDir=${denovoVCDir}/consensusHap
pafDir=${denovoVCDir}/paf
samDir=${denovoVCDir}/sam
vcfDir=${denovoVCDir}/vcf
vcfHg38Dir=${denovoVCDir}/vcfHg38
indelDir=${denovoVCDir}/indels
hetCalls=${denovoVCDir}/hetCalls
RefMHChaplotypes=${repoDir}/MHC_hapAssign/mhcRefHaps.txt
mkdir ${projectDir}
mkdir ${assemblyDir}
mkdir ${denovoVCDir}
mkdir ${consDir}
mkdir ${pafDir}
mkdir ${samDir}
mkdir ${vcfDir}
mkdir ${vcfHg38Dir}
mkdir ${indelDir}
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

conda activate amosPy27
conda env list
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
	 sh readQC.sh ${R1} ${R1A} ${R2} ${R2A} ${sampleAssembly} ${nThreads} ${primerFile} ${readCount} ${i} ${fastqDir} >> $readLog
	 cd ..
#########################################################################
	######## 1. Assign closest 2 reference MHC haplotypes ###########
	## ** Requires user to have HLA and C4 genotype files .txt **
        ## A) mHC_hapAssign(C4 genotypes, HLA genotypes) --> .txt
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

	####### 3. Error-correct and scaffold assembly ##########
	cd ./MHC_scaffold/
	sh refGuidedDeNovoAssembly_velvet_athena_scaffold_v3_target_haploid.sh $i ${R1A} ${R2A} ${fastqDir} ${sampleAssembly} ${insertSize} ${repoDir} ${Hap1} ${Hap2} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount} 
        
	if [ ${Hap1} != ${Hap2} ]
	then

		sh refGuidedDeNovoAssembly_velvet_athena_scaffold_v3_target_haploid.sh $i ${R1A} ${R2A} ${fastqDir} ${sampleAssembly} ${insertSize} ${repoDir} ${Hap2} ${Hap1} ${progBowtie2} ${progSamtools} ${progPicard} ${readCount}
	fi
	cd ..

	source activate py35
	conda env list

	###### 4. Scaffold and order assembly against refHap(s) #######
	cd ./MHC_order/
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
	sh run_denovo_varCall_samtools.sh ${i} ${Hap1} ${denovoVCDir}/consensusHap/${i}_vs_${Hap1}.fasta ${denovoVCDir} ${insertSize} ${sampleAssembly} ${progSamtools} 2> ${denovoVCDir}/indels/${i}_vs_${Hap1}_dn_indels.txt	
        
	if [ ${Hap1} != ${Hap2} ]
        then
                sh run_denovo_varCall_samtools.sh ${i} ${Hap2} ${denovoVCDir}/consensusHap/${i}_vs_${Hap2}.fasta ${denovoVCDir} ${insertSize} ${sampleAssembly} ${progSamtools} 2> ${denovoVCDir}/indels/${i}_vs_${Hap2}_dn_indels.txt
	fi
	cd ..


	######## 6. Merge scaffolded assembly vs hg38 var calls into 1 .vcf if heterozygous #############
	cd ./MHC_varHg38
	vcfExt=scaff_vs_hg38.vcf
	if [ ${Hap1} != ${Hap2} ] && [ ${Hap1} != 'Novel' ] && [ ${Hap1} != 'X' ] && [ ${Hap2} != 'Novel' ] && [ ${Hap2} != 'X' ] 
    then
		echo Merging VCF for sample $i
        V1="$(find ${vcfHg38Dir} -name ${i}*${Hap1}*${vcfExt})"
        V2="$(find ${vcfHg38Dir} -name ${i}*${Hap2}*${vcfExt})"
		echo $V1 $V2
        python newMergeVcf.py $V1 $V2 ${vcfHg38Dir}/${i}_${Hap1}_${Hap2}_mergedHet.vcf $Hap1 $Hap2
		echo Removing .gz and .csi file
		rm ${vcfHg38Dir}/*.gz
		rm ${vcfHg38Dir}/*.csi
		echo Done merging VCFs for sample $i
    fi
	cd ..
done
