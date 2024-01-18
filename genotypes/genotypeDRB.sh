#!/bin/bash
## Wrapper script to obtain the DRB genotypes needed for MHConstructor to run
## Usage: sh genotypeDRB.sh
## The fastq files location can be specified in ../control.txt
## The sample ID to run can be modified in ../testID.txt
## Output: a directory 'alleles' with all the T1K result and 'HLAgenotypes.csv' which processed
##         the result to a format that fits MHConstructor

source ../control.txt
repoDir=${binDir}/MHConstructor
progT1K=${repoDir}/tools/T1K/run-t1k
ids="$(cat ${repoDir}/testID.txt)"
set -e

for i in $ids; do
    echo $i
    R1="$(find ${targetCapture} -name ${i}*${R1ext})"
    R2="$(find ${targetCapture} -name ${i}*${R2ext})"
    if [[ -z "$R1" ]]; then
        continue
    fi

    ${progT1K} -1 ${R1} -2 ${R2} --preset hla -f ${repoDir}/tools/T1K/DRB1_IMGT_010324.fasta -o ${i} -t 8 --skipPostAnalysis --od ./alleles
    cd ./alleles
    find . ! -name '*_allele.tsv' -type f -exec rm -f {} +
    cd ../
done

python prepDRB.py ${repoDir}/genotypes/alleles