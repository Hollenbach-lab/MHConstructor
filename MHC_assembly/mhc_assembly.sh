#!/bin/bash

## Fastq processing- adapted from:
#########################################################
#  Reference-guided de novo assembly - SOAP
# ====================================================
# by Heidi Lischer, 2015/2016
#########################################################

## version 3
## 08/08/23

##To initialize conda env:
#conda activate amosPy27

#progPath=${8}
progPath=${15}

workPathFiles=${5}
readQCDir=${7}
# set work path ---------------------------
workPath=${7}/${1}_${4}_${6}

# log file
log=${workPath}/log_${1}_MHCassembly.txt

mkdir ${workPath}


# set variables #########################################
ref=${15}/refMHCh/${6}.fasta
refRed=${15}/refMHC/${6}

primerFile=${workPathFiles}/testPrimers.fasta
primerFileMP=${workPathFiles}/testPrimers.fasta

NThreads=8      # set the number of threads of every parallelizable step
maxReadLength=155
kmer=${12}         #define best K
declare -i expCov=${13}

# paired-end libraries -------------------
name=$1           # set name of your species
lib=(150)      # set insertion libraries
insSize=(${11})
insLow=(0)     # lower bound of insertion size
libsd=(100)       # sd of insertion size

# list of files with forward reads according to lib array
reads1=(${3})
# list of files with rewerse reads according to lib array
reads2=(${4})
# short names of libraries
shortNames=(${1})

# Programs -------------------------------- 
progIdba=idba
progFastQC=fastqc ##conda
progTrimmomatic=trimmomatic ##conda
#progSamtools=${progPath}/samtools-1.16.1/samtools
progSamtools=${9}
progVcfutils=vcfutils.pl
#progBcftools=${progPath}/bcftools-1.16.1/bcftools
progBcftools=bcftools
progBamtools=bamtools ##conda
progBedtools=bedtools ##conda
progPicard=${10}
#progBowtie2=/home/kwade/bin/bowtie2-2.4.2-sra-linux-x86_64/bowtie2
progBowtie2=${16}
progSeqtk=seqtk ##conda
progNucmer=nucmer ## installed through 'mummer'
progGatk=gatk ##conda
progSoapdenovo2=soapdenovo2 ##conda
progRemovShortSeq=${15}/MHC_assembly/RemoveShortSeq.jar
progGetBlocks=${15}/MHC_assembly/GetBlocks.jar
progFastaToAmos=${15}/MHC_assembly/FastaToAmos.jar
progWriteSoapConfig=${15}/MHC_assembly/WriteSoapConfig.jar
progFastaStats=${15}/MHC_assembly/FastaStats.jar
progSplitSeqLowCov=${15}/MHC_assembly/SplitSeqLowCov.jar


# run reference-guided de novo MHC assembly ##########################################

## Initialize processed fastq read variables ####
  cd $workPath
  read1TrimPair=()
  read1TrimUnPair=()
  read2TrimPair=()
  read2TrimUnPair=()
  for i in ${!lib[*]}  #for all indexes in the array
  do
    read1TrimPair[i]=${readQCDir}/${shortNames[i]}_R1_trimPair.fastq
    read1TrimUnPair[i]=${readQCDir}/${shortNames[i]}_R1_trimUnPair.fastq
    read2TrimPair[i]=${readQCDir}/${shortNames[i]}_R2_trimPair.fastq
    read2TrimUnPair[i]=${readQCDir}/${shortNames[i]}_R2_trimUnPair.fastq
  done
  readTrimUnPair=${readQCDir}/${name}_trimUnpair_mod.fastq
  cat ${read1TrimUnPair[*]} ${read2TrimUnPair[*]} > ${readTrimUnPair}
  libUnpair=Unpair

# 2. Step: map reads against reference
#          and define blocks and superblocks
#######################################################
  # prepare Reference ---------
  echo "prepare reference..."
  echo "prepare reference..." >> $log

  #remove scaffolds shorter than 10 kb
  java -jar ${progRemovShortSeq} -i $ref -o ${refRed}.fa -length 1000

  #create index files
  ### *** Need to replace when picard is integrated ** #####
  ##${progSamtools} faidx ${refRed}.fa
  ##${progPicard} CreateSequenceDictionary R=${refRed}.fa O=${refRed}.dict

  # map reads against reference ----------
  echo "run reference mapping..."
  echo "run reference mapping..." >> $log

  #cd ${workPath}

  # index reference file
  echo "run reference mapping..."
  
  ${progBowtie2}-build ${refRed}.fa ${refRed}

  mappedAll=()
  unmapped=()
  mapped=()
  mappedFiltered=()
  bowtieFailPair=()
  bowtieFailUnPair1=()
  bowtieFailUnPair2=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do
    mappedAll[i]=${workPath}/${shortNames[i]}_all.sorted.bam
    unmapped[i]=${workPath}/${shortNames[i]}_unmapped.sorted.bam
    mapped[i]=${workPath}/${shortNames[i]}.sorted.bam
    mappedFiltered[i]=${workPath}/${shortNames[i]}.sorted.filtered.bam
    bowtieFailPair[i]=${workPath}/${shortNames[i]}_failPair.fastq
    bowtieFailUnPair1[i]=${workPath}/${shortNames[i]}_failUnPairR1.fastq
    bowtieFailUnPair2[i]=${workPath}/${shortNames[i]}_failUnPairR2.fastq
    bowtieFailUnpairMerged[i]=${workPath}/${shortNames[i]}_failUnP.fastq    
    (
      ${progBowtie2} --fast-local --threads 8 -q --phred33 -I ${insLow} -X ${insSize} -x ${refRed} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} -U ${read1TrimUnPair[i]},${read2TrimUnPair[i]} | ${progSamtools} view -bS | ${progSamtools} sort -T ${shortNames[i]} -o ${mappedAll[i]}
      ${progSamtools} index ${mappedAll[i]}

      #filter unmapped reads
      ${progSamtools} view -b -F 4 ${mappedAll[i]} > ${mapped[i]}
      ${progSamtools} index ${mapped[i]}

      #get unmapped reads
      ${progSamtools} view -b -f 4 ${mappedAll[i]} > ${unmapped[i]}
      ${progSamtools} view -b -f 9 ${unmapped[i]} > ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progPicard} SamToFastq I=${unmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${bowtieFailPair[i]%.fastq}.1.fastq F2=${bowtieFailPair[i]%.fastq}.2.fastq
      rm ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progSamtools} view -b -F 8 -f 64 ${unmapped[i]} | ${progBamtools} convert -format fastq -out ${bowtieFailUnPair1[i]}
      ${progSamtools} view -b -F 8 -f 128 ${unmapped[i]} | ${progBamtools} convert -format fastq -out ${bowtieFailUnPair2[i]}
      ${progBamtools} stats -in ${mappedAll[i]} >> $log
      echo "--> ${mappedAll[i]}" >> $log
      cat ${bowtieFailUnPair1[i]} ${bowtieFailUnPair2[i]} > ${bowtieFailUnpairMerged[i]}
      #filter for mapping quality >=10
      ${progSamtools} view -b -q 10 ${mapped[i]} > ${mappedFiltered[i]}
      ${progBamtools} stats -in ${mappedFiltered[i]} >> $log
      echo "--> ${mappedFiltered[i]}" >> $log

      #check insertion size
      ${progPicard} CollectInsertSizeMetrics R=${refRed}.fa I=${mapped[i]} O=${mapped[i]%.bam}_insertSize.txt H=${mapped[i]%.bam}_insertSizeHist.pdf
      ) &
      let count+=1
      [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
  
  #merge alignment files
  mappedMerged=${workPath}/${name}.sorted_mapped
  #rm ${mappedMerged}.bam
  ${progSamtools} merge -f ${mappedMerged}.bam ${mapped[*]}
  #${progSamtools} merge -f ${mappedMerged}.bam ${mapped[*]} ${mateMapped[*]}
  # get blocks and superblocks ----------
  echo "get blocks and superblocks..."
  echo "get blocks and superblocks..." >> $log

  #get coverage along genome
  covFile=${workPath}/${name}_coverage.txt
  ${progBedtools} genomecov -ibam ${mappedMerged}.bam -bga > ${covFile}
  #${progBedtools} genomecov -ibam ${mappedMerged}.bam -bga > ${covFile}
  #only proparly paired reads
  ${progSamtools} view -bf 0x2 ${mappedMerged}.bam | ${progSamtools} sort - -n | ${progBedtools} bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | ${progBedtools} genomecov -i - -bga -g ${refRed}.fa.fai > ${covFile%.txt}Paired.txt

  #get blocks with a minimal coverage of 10 (paired-end) reads and create superblocks of at least 12000bp length and a minimal overlap of 300bp (max. overlap = 3*300bp)
  blocks=${workPath}/blocks.txt
  superblocks=${workPath}/superblocks.txt
  java -jar ${progGetBlocks} -i ${covFile} -paired ${covFile%.txt}Paired.txt -o ${blocks} -oSuper ${superblocks} -mCov 10 -sLength 12000 -sOverlap 300 -maxLength 10000000


# 3. Step: do deNovo assembly within superblocks
#######################################################
  echo "deNovo assembly within superblocks..."
  echo "deNovo assembly within superblocks..." >> $log
  #source /home/kwade/.conda/envs/assembleMHC/bin/

  cd ${workPath}
  velvetRes=${workPath}/velvetResults
  mkdir ${velvetRes}

  #split superbolcks.txt into ${NThreads} files to run it in parallel
  size=$(($(wc -l < ${superblocks})/$((NThreads))+1))
  split -l $size -d ${superblocks} ${superblocks%.txt}_

  count=0
  for (( j=0; j<$((NThreads)); j++ ))
  do
    file=${superblocks%.txt}_0${j}
    fileout=${superblocks%.txt}_0${j}_run.sh
    logout=${superblocks%.txt}_0${j}.log

    array=(${file//_/ })
    number=${array[${#array[*]}-1]}
    if [ ${number} != "00" ]
    then
      number=`echo $number|sed 's/^0*//'`
    fi

    if [[ $number =~ ^[0-9]+$ ]]
    then
      blockNb=$(($number*$size+1))
    else
      blockNb=1
    fi

    printf "#"'!'"/bin/bash\n"                 > ${fileout}
    printf "\n"                                >> ${fileout}
    printf "mkdir ${file}_temp\n"              >> ${fileout}
    printf "cd ${file}_temp\n"                 >> ${fileout}
    printf "\n" >> ${fileout}

    printf "blockNb=$blockNb\n" >> ${fileout}
    printf "start=\`date +%%s\`\n" >> ${fileout}
    printf "for block in \$(cat ${file})\n" >> ${fileout}
    printf "do\n" >> ${fileout}
    printf "  echo \$blockNb >> $logout\n" >> ${fileout}
    printf "\n" >> ${fileout}
    printf "  #extract sequence names within specified region\n" >> ${fileout}
    #block="$(sed -n "${blockNb}p" $file)"

    echo "Added #############"
    #echo "$(\$block)"
    seqNames=()
    seqBam=()
    subSeq1=()
    subSeq2=()
    blockSize=$(less $file| wc -l)

    for i in ${!lib[*]}
    do
      seqNames[i]=sequences_${shortNames[i]}
      seqBam[i]=sequences_${shortNames[i]}.bam
      subSeq1[i]=subseq_${shortNames[i]}_R1.fastq
      subSeq2[i]=subseq_${shortNames[i]}_R2.fastq

      #1-extract reads where both pairs mapped
      printf "  ${progSamtools} view -b ${mapped[i]} \$block | ${progSamtools} sort -T ${shortNames[i]} -n -o ${seqBam[i]}\n" >> ${fileout}
      ## Edit from original: using samtools fastq instead bcftools bamtofastq
      printf " ${progSamtools} fastq -1 1_${subSeq1[i]} -2 1_${subSeq2[i]} ${seqBam[i]}\n" >> ${fileout}
     
      #extract paired reads with one pair unmapped
      #2- First in pair didn't map, other did
      ### 10/19/22- Correction: from samtools doc: 72- mapped read is first in pair, other didn't map
      printf "  ${progSamtools} view  -b -f 72 ${seqBam[i]} | ${progBamtools} convert -format fastq | paste - - - - | sort -k 1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 2_${subSeq1[i]}\n" >> ${fileout}
      ##  10/19/22-Switched to '2.2'
      printf "  ${progSamtools} view  -f 72 ${seqBam[i]} | cut -f 1 | rev | cut -c2- | rev |awk \'{print \$0\"2.2\"}\' > ${seqNames[i]}_R2.txt\n" >> ${fileout}
      ### 10/19/22- Switched to 1.1
      printf "  ${progSamtools} view  -f 136 ${seqBam[i]} | cut -f 1 | rev | cut -c2- | rev |awk \'{print \$0\"1.1\"}\' > ${seqNames[i]}_R1.txt\n" >> ${fileout}
      ### 10/19/22- Correction: from samtools doc: 136- mapped read is second in pair, other didn't map
      #3- extract paired reads with second in pair unmapped, other did
      printf "  ${progSamtools} view  -b -f 136 ${seqBam[i]} | ${progBamtools} convert -format fastq | paste - - - - | sort -k 1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 3_${subSeq2[i]}\n" >> ${fileout}

      printf "  ${progSeqtk} subseq ${bowtieFailUnPair2[i]} ${seqNames[i]}_R2.txt | paste - - - - | sort -k 1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 2_${subSeq2[i]}\n" >> ${fileout}
      printf "  ${progSeqtk} subseq ${bowtieFailUnPair1[i]} ${seqNames[i]}_R1.txt | paste - - - - | sort -k 1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 3_${subSeq1[i]}\n" >> ${fileout}

      printf "  cat 1_${subSeq1[i]} 2_${subSeq1[i]} 3_${subSeq1[i]} > ${subSeq1[i]}\n" >> ${fileout}
      printf "  cat 1_${subSeq2[i]} 2_${subSeq2[i]} 3_${subSeq2[i]} > ${subSeq2[i]}\n" >> ${fileout}
      printf "\n" >> ${fileout}
    done


    printf "  #velvet --------\n" >> ${fileout}
    subSeqFasta=()
    for i in ${!lib[*]}
    do
      subSeqFasta[i]=subseq_${lib[i]}.fa
      printf "  fq2fa --merge ${subSeq1[i]} ${subSeq2[i]} ${subSeqFasta[i]}\n" >> ${fileout}
      ### 8/23-Unpaired reads are already included in the R1, R2 subseq files
      printf "  mv ${subSeqFasta[i]} subseq.fa\n" >> ${fileout}
    done

    subSeqFasta[${!subSeqFasta[*]}]=subseq_UnPair.fa
    printf "  velveth velvetResults ${kmer} -shortPaired subseq.fa\n" >> ${fileout}
    
    printf "  velvetg velvetResults -exp_cov ${expCov} -ins_length ${insSize}\n" >> ${fileout}
    printf "  if [ ! -f velvetResults/contigs.fa ]\n" >> ${fileout}
    printf "  then\n" >> ${fileout}
    printf "    echo \"\$blockNb velvet failed\" >> $logout\n" >> ${fileout}
    printf "  fi\n" >> ${fileout}
    printf "  cp velvetResults/contigs.fa ${velvetRes}/sblock\${blockNb}-contigs.fa\n" >> ${fileout}
    printf "\n" >> ${fileout}

    printf "  rm -rf ${file}_temp/*\n" >> ${fileout}
    printf "  ((blockNb++))\n" >> ${fileout}
    printf "done\n" >> ${fileout}
    printf "\n" >> ${fileout}
    printf "rm -rf ${file}_temp\n" >> ${fileout}
    printf "end=\`date +%%s\`\n" >> ${fileout}
    printf "echo \$((end-start))\n" >> ${fileout}
    printf "\n" >> ${fileout}
   chmod +x ${fileout}
    ${fileout} &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait

# deNovo assembly of unassembled reads ----------
  echo "deNovo assembly of unassembled reads..."
  echo "deNovo assembly of unassembled reads..." >> $log

  unassFolder=${workPath}/Unassembled
  mkdir ${unassFolder}
  cd ${unassFolder}

  subSeqFasta=()

  for i in ${!lib[*]}
  do
    subSeqFasta[i]=${shortNames[i]}_failPair.fa
    fq2fa --merge ${bowtieFailPair[i]%.fastq}.1.fastq ${bowtieFailPair[i]%.fastq}.2.fastq ${subSeqFasta[i]}
  done
  subSeqFasta[${#subSeqFasta[*]}]=${name}_failUnp.fa
  fq2fa ${bowtieFailUnpairMerged[*]} ${subSeqFasta[${#subSeqFasta[*]}-1]}
  cat ${subSeqFasta[*]} > subseq.fa

  ## Replacing with hard code expected coverage for now 8/23/23 ##
  ##expCovLow=`expr $expCov / 2`
  velveth velvetResults ${kmer} -shortPaired subseq.fa
  #velvetg velvetResults -exp_cov $expCovLow -ins_length ${insSize}
  velvetg velvetResults -exp_cov 30 -ins_length ${insSize}



# 4. Step: get non-redundant supercontigs
####################################################### 
  echo "get supercontigs..."
  echo "get supercontigs..." >> $log
  cd ${workPath}
  
 

  #merge deNovo assembled superblocks into one FASTA file
  velvetContigs=${velvetRes}/velvetContigs.fa
  cat ${velvetRes}/*-contigs.fa > ${velvetContigs}

  #merge Unassembled files
  orgUnass=${unassFolder}/velvetResults/contigs.fa
  unass500=${unassFolder}/Unassembled_velvet_500.fa
  java -jar ${progRemovShortSeq} -i ${orgUnass} -o ${unass500} -length 500  >> $log
  
  ### Map to hg38 genome and alt MHC hap (if heterozygous) and remove off-target mapping, unaligned contigs ####
  mkdir ${unassFolder}/offTarget
  cd ${15}/MHC_assembly
  sh checkAndRemoveOfftarget.sh ${shortNames} hg38Patch11 ${unass500} ${unassFolder} ${14} ${6}
  unass500OT=${unassFolder}/Unassembled_velvet_OT.fa
  cd ${workPath}
  ###################

#merge all files
  superblockSeq=${workPath}/deNovo_Superblocks.fa
  cat ${velvetContigs} ${unass500OT} > ${superblockSeq}

  #remove short seq (<200)
  echo "remove seq < 200 in ${superblockSeq}" >> $log


  superblockSeq200=${superblockSeq%.fa}_100.fa
  java -jar ${progRemovShortSeq} -i ${superblockSeq} -o ${superblockSeq200} -length 100 -u  >> $log


  #remove redundency with AMOScmp
  amosFolder=${workPath}/AMOScmp
  mkdir ${amosFolder}
  cd ${amosFolder}

  #assemble all assembled superblocks with AMOScmp to supercontigs (with the help of reference)
  #changed parameters in AMOScmp: (casm-layout -t 1000 (maximum ignorable trim length), make-consensus -o 10 (minimum overlap base))
  superblockSeqAmos=${superblockSeq200%.fa}_Amos.afg
  java -jar ${progFastaToAmos} -i ${superblockSeq200} -o ${superblockSeqAmos}
  supercontigs=Amos_supercontigs
  amosSupercontigs=${amosFolder}/${supercontigs}.fasta

  echo "run AMPScmp..." >> $log
  
  rm -r ./Amos_supercontigs.bnk
 
  ################################

  # running AMPScmp step by step and use multithread nucmer to spead it up
  ## Building AMOS bank
  echo "  build AMPS bank..." >> $log
  bank-transact -c -z -b ${supercontigs}.bnk -m ${superblockSeqAmos}

  bank-unlock Amos_supercontigs.bnk
  ## Collecting clear range sequences
  echo "  clear range sequences..." >> $log
  dumpreads ${supercontigs}.bnk > ${supercontigs}.seq

  ## Running nucmer
  echo "  run nucmer..." >> $log
  ${progNucmer} --maxmatch --prefix=${supercontigs} ${refRed}.fa ${supercontigs}.seq
  
  ##ADDED 9/26/22
  bank-unlock Amos_supercontigs.bnk

  ## Running layout
  echo "  run layout..." >> $log
  casm-layout -t 1000 -U ${supercontigs}.layout -C ${supercontigs}.conflict -b ${supercontigs}.bnk ${supercontigs}.delta


  ## Running consensus
  echo "  run consensus..." >> $log
  make-consensus -o 10 -B -b ${supercontigs}.bnk
 

  ## Outputting contigs
  echo "  output contigs..." >> $log
  bank2contig ${supercontigs}.bnk > ${supercontigs}.contig

  ## Outputting fasta
  echo "  output fasta..." >> $log
  bank2fasta -b ${supercontigs}.bnk > ${supercontigs}.fasta


# 5. Step: map reads on supercontigs
#          and de novo assemble unmapped reads
####################################################### 
  echo "map reads on supercontigs and correct them..."
  echo "map reads on supercontigs and correct them..." >> $log

  #make seqnames unique
  amosSupercontigsUnique=${amosSupercontigs%.fasta}_unique.fa


  java -jar ${progRemovShortSeq} -i ${amosSupercontigs} -o ${amosSupercontigsUnique} -length 1 -u

  #get statistics
  echo ${amosSupercontigsUnique} >> $log
  java -jar ${progFastaStats} -i ${amosSupercontigsUnique} -min 100 >> $log

