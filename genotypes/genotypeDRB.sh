#!/bin/bash
## Wrapper script to obtain the DRB genotypes needed for MHConstructor to run
## Usage: sh genotypeDRB.sh
## The fastq files location can be specified in ../control.txt
## The sample ID to run can be modified in ../testID.txt

source ../control.txt
repoDir=${binDir}/MHConstructor
progT1K=${repoDir}/tools/T1K/run-t1k
ids="$(cat ${repoDir}/testID.txt)"
set -e

for i in $ids; do
    echo $i
    ${progT1K} -1 ${targetCapture}/${i}_${R1ext} -2 ${targetCapture}/${i}_${R2ext} --preset hla -f ${repoDir}/tools/T1K/DRB1_IMGT_010324.fasta -o ${i} -t 8 --skipPostAnalysis --od ./alleles
    cd ./alleles
    find . ! -name '*_allele.tsv' -type f -exec rm -f {} +
    cd ../
done

python prepDRB.py ${repoDir}/genotypes/alleles