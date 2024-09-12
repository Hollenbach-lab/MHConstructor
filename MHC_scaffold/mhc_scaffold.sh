

#!/bin/bash

## Adapted from:
#########################################################
#  Reference-guided de novo assembly - SOAP
# ====================================================
# by Heidi Lischer, 2015/2016
# https://bitbucket.org/HeidiLischer/refguideddenovoassembly_pipelines/src/master/
#########################################################
## Edited: 9/11/24
## Adapted by: wadekj and rsuseno

workPathFiles=${4}

# set work path ---------------------------
workPath=${5}

# log file
log=${workPath}/log_${1}.txt


# set variables #########################################
HapA=${7}/MHC_generate/refMHChaps/${8}
HapB=${7}/MHC_generate/refMHChaps/${9}



NThreads=8      # set the number of threads of every parallelizable step
maxReadLength=155
kmer=${13}	

# paired-end libraries -------------------
name=${1}           # set name of your species
lib=(150)      # set insertion libraries
insLow=(0)     # lower bound of insertion size
libsd=(400)       # sd of insertion size

# list of files with forward reads according to lib array
reads1=(${2})
# list of files with rewerse reads according to lib array
reads2=(${3})
# short names of libraries
shortNames=${1}



# Programs -------------------------------- 
progPath=${7}
progIdba=idba
progFastQC=fastqc ##conda
progTrimmomatic=trimmomatic ##conda
progSamtools=${11}
progVcfutils=vcfutils.pl
progBcftools=${progPath}/tools/bcftools-1.16/bcftools
progBamtools=bamtools ##conda
progBedtools=bedtools ##conda
progPicard=${12}
progBowtie2=${10}
progSeqtk=seqtk ##conda
progNucmer=nucmer ## installed through 'mummer'
progGatk=gatk
progFusion=${progPath}/tools/SOAPdenovo2/
progRemovShortSeq=${progPath}/MHC_assembly/RemoveShortSeq.jar
progGetBlocks=${progPath}/MHC_assembly/GetBlocks.jar
progFastaToAmos=${progPath}/MHC_assembly/FastaToAmos.jar
progWriteSoapConfig=${progPath}/MHC_assembly/WriteSoapConfig.jar
progFastaStats=${progPath}/MHC_assembly/FastaStats.jar
progSplitSeqLowCov=${progPath}/MHC_assembly/SplitSeqLowCov.jar

### Initialize read file variables #####
  read1TrimPair=()
  read1TrimUnPair=()
  read2TrimPair=()
  read2TrimUnPair=()
  for i in ${!lib[*]}  #for all indexes in the array
  do
    read1TrimPair[i]=${workPath}/${shortNames[i]}_R1_trimPair.fastq
    read1TrimUnPair[i]=${workPath}/${shortNames[i]}_R1_trimUnPair.fastq
    read2TrimPair[i]=${workPath}/${shortNames[i]}_R2_trimPair.fastq
    read2TrimUnPair[i]=${workPath}/${shortNames[i]}_R2_trimUnPair.fastq
  done
###################################
# 5. Step: map reads on supercontigs
#          and de novo assemble unmapped reads
####################################################### 
  echo "map reads on supercontigs and correct them..."
  echo "map reads on supercontigs and correct them..." >> $log

  amosFolder=${5}/${1}_${13}_${8}/AMOScmp
  mkdir $amosFolder
  cd $amosFolder
  amosSupercontigsUnique=${amosFolder}/Amos_supercontigs_unique.fa

#prepare reference
  ${progBowtie2}-build ${amosSupercontigsUnique} ${amosSupercontigsUnique%.fa}

  supercontMappedAll=()
  supercontUnmapped=()
  supercontFailPair=()
  supercontFailUnpair=()
  supercontMappedFiltered=()
  count=0

  read1=${workPathFiles}/${reads1}
  read2=${workPathFiles}/${reads2}
  for i in ${!lib[*]}  #for all indexes in the array
  do
    supercontMappedAll[i]=${amosFolder}/${shortNames[i]}_all.sorted.bam
    supercontUnmapped[i]=${amosFolder}/${shortNames[i]}_unmapped.sorted.bam
    supercontFailPair[i]=${amosFolder}/${shortNames[i]}_failPair.fastq
    supercontFailUnpair[i]=${amosFolder}/${shortNames[i]}_failUnp.fastq
    supercontMappedFiltered[i]=${amosFolder}/${shortNames[i]}.filtered.sorted.bam
    (
       ${progBowtie2} --sensitive -p 8 -q --phred33 -I 0 -X ${6} -x ${amosSupercontigsUnique%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${shortNames[i]} -o ${supercontMappedAll[i]}
      ${progSamtools} index ${supercontMappedAll[i]}

      #get unmapped reads
      ${progSamtools} view -b -f 4 ${supercontMappedAll[i]} > ${supercontUnmapped[i]}
      ${progSamtools} view -b -f 9 ${supercontUnmapped[i]} > ${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progPicard} SamToFastq I=${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${supercontFailPair[i]%.fastq}.1.fastq F2=${supercontFailPair[i]%.fastq}.2.fastq
      ${progSamtools} view -b -F 8 ${supercontUnmapped[i]} | ${progSamtools} bam2fq > ${supercontFailUnpair[i]}
            ${progBamtools} stats -in ${supercontMappedAll[i]} >> $log
      echo "--> ${supercontMappedAll[i]}" >> $log

      #filter for mapping quality >=10
      ${progSamtools} view -b -F 4 -q 1 ${supercontMappedAll[i]} > ${supercontMappedFiltered[i]}
      ${progBamtools} stats -in ${supercontMappedFiltered[i]} >> $log
      echo "--> ${supercontMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait



  # 6. Step: Correct supercontigs
##########################################################
  echo "merge contigs..."
  echo "merge contigs..." >> $log
  mergedFolder=${workPath}/${1}_${13}_${8}/merged_corr
  mkdir $mergedFolder
  cd $mergedFolder

  ## 3/1/23: replace
  ##merged=${mergedFolder}/${name}_supercontSeq_Unass.fa
  ##cat ${amosSupercontigsUnique} ${supercontSeqUnass} > $merged
  merged=${amosSupercontigsUnique}

  echo "${merged}" >> $log
  java -jar ${progFastaStats} -i ${merged} -min 500 >> $log

  ${progBowtie2}-build ${merged} ${merged%.fa}
 
  mergedMappedAll=()
  mergedMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do
    read1TrimPair[i]=${workPath}/${shortNames[i]}_R1_trimPair.fastq
    read2TrimPair[i]=${workPath}/${shortNames[i]}_R2_trimPair.fastq
    read1original[i]=${originalR1[i]}
    read2original[i]=${originalR2[i]}
    mergedMappedAll[i]=${mergedFolder}/${shortNames[i]}_all.sorted.bam
    mergedMappedFiltered[i]=${mergedFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      ${progBowtie2} --sensitive -p 8 -q --phred33 -I 0 -X ${6} -x ${merged%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${shortNames[i]} -o ${mergedMappedAll[i]}
      ${progSamtools} index ${mergedMappedAll[i]}

      ${progBamtools} stats -in ${mergedMappedAll[i]} >> $log
      echo "--> ${mergedMappedAll[i]}" >> $log

      #filter for mapping quality >=10
      ${progSamtools} view -b -F 4 -q 1 ${mergedMappedAll[i]} > ${mergedMappedFiltered[i]}

            ${progBamtools} stats -in ${mergedMappedFiltered[i]} >> $log
      echo "--> ${mergedMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait


  # error correction ----------
  # add RG header of second file (lost while merging)
  for i in ${!lib[*]}  #for all indexes in the array
  do
    echo -e "@RG\tID:${shortNames[i]}.filtered.sorted\tPL:illumina\tPU:${lib[i]}\tLB:${lib[i]}\tSM:${shortNames[i]}" >> rg
  done
  ${progSamtools} view -H ${mergedMappedFiltered[*]} | cat - rg > header
  
 
  #realign reads
  mergedMappedMerged=${mergedFolder}/${name}.filtered_RG.sorted.bam
  rm $mergedMappedMerged
  ${progSamtools} merge -r -h header ${mergedMappedMerged} ${mergedMappedFiltered[*]}
  rm rg header
  ${progSamtools} index ${mergedMappedMerged}
  ${progSamtools} faidx ${merged}
  rm ${merged%.fa}.dict
  ${progPicard} CreateSequenceDictionary R=${merged} O=${merged%.fa}.dict
  
  ##### ****** MAJOR ANALYSIS DECISION ******** ############
  ## 9/26 Phased out of GATK v4- HaplotypeCaller is preferable to use now
  ## Replacing with haplotypecaller
  ## *********************************************

  mergedMappedMergedReal=${mergedMappedMerged%.bam}_realigned.bam
  
  ${progGatk} HaplotypeCaller -R ${merged} -I ${mergedMappedMerged} -O ${name}supercontigs.vcf.gz -bamout ${mergedMappedMergedReal}  
######################################3 
 
  #get alternative seq
  mergedCorr=${merged%.fa}_corr.fa

  #### New 3/2/23- updated from vcfutil vcf2fq to bcftools consensus ####
  
  ## Call variants
  ${progBcftools} mpileup -Ou -d 100000 -f ${merged} ${mergedMappedMergedReal} | ${progBcftools} call -mv -Oz -o ${merged%.fa}_corr.vcf.gz
  
  ## Index calls
  ${progBcftools} index ${merged%.fa}_corr.vcf.gz

  ## Normalize indels
  ${progBcftools} norm -f ${merged} ${merged%.fa}_corr.vcf.gz -Ob -o ${merged%.fa}_corr.norm.bcf
  
  ## filter adjacent indels within 5bp
  ${progBcftools} filter --IndelGap 5 ${merged%.fa}_corr.norm.bcf -Ob -o ${merged%.fa}_corr.norm.flt-indels.bcf
  
  ## apply variants to create consensus seq
  cat ${merged} | ${progBcftools} consensus ${merged%.fa}_corr.vcf.gz  > ${mergedCorr}
  
  mergedCorrWN=${mergedCorr}

  #get statistics
  echo ${mergedCorrWN} >> $log
  java -jar ${progFastaStats} -i ${mergedCorrWN} -min 500 >> $log


  # split sequences at places with no coverage ----------
  ${progBowtie2}-build ${mergedCorrWN} ${mergedCorrWN%.fa}

  mergedCorrMappedAll=()
  mergedCorrMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do
      read1TrimPair[i]=${workPath}/${shortNames[i]}_R1_trimPair.fastq
      read2TrimPair[i]=${workPath}/${shortNames[i]}_R2_trimPair.fastq
      read1original[i]=${originalR1[i]}
      read2original[i]=${originalR2[i]}
      mergedCorrMappedAll[i]=${mergedFolder}/${shortNames[i]}_corrWN_all.sorted.bam
      mergedCorrMappedFiltered[i]=${mergedFolder}/${shortNames[i]}_corrWN.filtered.sorted.bam
      ${progBowtie2} --sensitive -p 8 -q --phred33 -I 0 -X ${6} -x ${mergedCorrWN%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${shortNames[i]} -o ${mergedCorrMappedAll[i]}
      ${progSamtools} index ${mergedCorrMappedAll[i]}

      ${progBamtools} stats -in ${mergedCorrMappedAll[i]} >> $log
      echo "--> ${mergedCorrMappedAll[i]}" >> $log

      #filter for mapping quality >=10
     
      ${progSamtools} view -b -F 4 -q 1 ${mergedCorrMappedAll[i]} > ${mergedCorrMappedFiltered[i]}

      ${progBamtools} stats -in ${mergedCorrMappedFiltered[i]} >> $log
      echo "--> ${mergedCorrMappedFiltered[i]}" >> $log
    
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait

 

  mergedCorrMappedFilteredMerged=${mergedFolder}/${name}_corrWN.filtered.sorted.bam
  
  ## 11/1
  ${progBedtools} genomecov -ibam ${mergedCorrMappedFilteredMerged} -bga > ${mergedCorrWN%.fa}_filteredCov.txt
  #only proparly paired reads
  ${progSamtools} faidx $mergedCorrWN
  ${progSamtools} view -bf 0x2 ${mergedCorrMappedFilteredMerged} | ${progSamtools} sort - -n | ${progBedtools} bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | ${progBedtools} genomecov -i - -bga -g ${mergedCorrWN}.fai  > ${mergedCorrWN%.fa}_filteredPairedCov.txt

  java -jar ${progSplitSeqLowCov} -i ${mergedCorrWN%.fa}_filteredCov.txt -paired ${mergedCorrWN%.fa}_filteredPairedCov.txt -o ${mergedCorrWN%.fa}_filteredNotCov.txt -mCov 1 -fasta ${mergedCorrWN} -fastaOut ${mergedCorrWN%.fa}_splitFiltered.fa >> $log
  echo ${mergedCorrWN%.fa}_splitFiltered.fa >> $log
  java -jar ${progFastaStats} -i ${mergedCorrWN%.fa}_splitFiltered.fa -min 100 >> $log


  # 7. Step: scaffolding and gap closing
#######################################################
  echo "scaffolding..."
  echo "scaffolding..." >> $log

  scafFolder=${mergedFolder}/scaffold_gapClosed
  mkdir $scafFolder
  cd $scafFolder

  for i in ${!lib[*]}
  do
    read1TrimPair[i]=${workPath}/${shortNames[i]}_R1_trimPair.fastq
    read2TrimPair[i]=${workPath}/${shortNames[i]}_R2_trimPair.fastq
    read1original[i]=${originalR1[i]}
    read2original[i]=${originalR2[i]}
    if [ $i == 0 ]
    then
      libList=${lib[i]}
      forwardReads=${read1TrimPair[i]}
      reverseReads=${read2TrimPair[i]}
     else
      libList=${libList},${lib[i]}
      forwardReads=${forwardReads},${read1TrimPair[i]}
      reverseReads=${reverseReads},${read2TrimPair[i]}
    fi
  done


  #write config file
  soapConf=${scafFolder}/soap.config
  java -jar ${progWriteSoapConfig} -insLength ${6} -r1 ${forwardReads} -r2 ${reverseReads} -max ${maxReadLength} -ru 2 -rank 1 -o ${soapConf}
  scafFile=${name}_${kmer}


  ${progFusion}SOAPdenovo-fusion -D -c ${mergedCorrWN%.fa}_splitFiltered.fa -K ${kmer} -g ${scafFile} -p ${NThreads}
  #SOAPdenovo-127mer map -s ${soapConf} -g ${scafFile} -p ${NThreads}
  #SOAPdenovo-127mer scaff -g ${scafFile} -p ${NThreads} -F
  
  #SOAPdenovo-127mer pregraph -s ${soapConf} -o ${scafFile} -R -K ${kmer} -p ${NThreads} 
  #SOAPdenovo-127mer contig -g ${scafFile} -p ${NThreads}
  SOAPdenovo-63mer map -s ${soapConf} -g ${scafFile} -p ${NThreads}
  SOAPdenovo-63mer scaff -g ${scafFile} -p ${NThreads}

  GapCloser -a ${scafFile}.scafSeq -b soap.config -o ${scafFile}_final.fasta -l 155 -t ${NThreads}

  #remove scaffolds < 100 bp ----------
  scafSeq=${scafFolder}/${name}_scafSeq.fa
  
  #echo ${scafFile}.scafSeq >> $log
  java -jar ${progRemovShortSeq} -i ${scafFile}_final.fasta -o ${scafFile}_final_100.fasta -length 100 >> $log

  #get statistics
  echo ${scafSeq} >> $log
  java -jar ${progFastaStats} -i ${scafFile}_final.fasta -min 100 >> $log
  java -jar ${progFastaStats} -i ${scafFile}_final.fasta -min 500 >> $log

  #map reads against scaffolds
  ${progBowtie2}-build ${scafFile}_final_100.fasta ${scafFile}_final_100
