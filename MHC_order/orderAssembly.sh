#! /bin/bash
## Scaffold assembly to guide BMH sequence using RagTag
## Edited: 9/11/24
## Authors: wadekj and susenor

# $1 - refHapID
# $2 - assembly fasta
# $3 - AssemblyDir/merged/scaffolds
# $4 - sampleID
# $5 - projectDir
ragtag.py scaffold ../refMHC/${2}.fasta ${1} -r -o ${3}/RT_${2}
ragtag.py patch ${3}/RT_${2}/ragtag.scaffold.fasta ../refMHC/${2}.fasta -o ${3}/RT_${2}

## Extract only the continuous scaffolded sequence consisting of placed contigs
head -n 2 ${3}/RT_${2} > ${5}/consensusHap/${4}_${2}_RT_onlyPlaced.fasta

## Extract the unplaced contigs to separate file
tail -n +3 ${3}/RT_${2}/ragtag.scaffold.fasta > ${5}/consensusHap/unplacedContigs/${4}_${2}_RT_unPlaced.fasta
