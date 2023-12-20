#! /bin/bash

# $1 - refHapID
# $2 - assembly fasta

#ragtag.py correct ../refHaps/${Hap1}
ragtag.py scaffold ../refMHChaps/${2}.fasta ${1} -r -o ${3}/RT_${2}
ragtag.py patch ${3}/RT_${2}/ragtag.scaffold.fasta ../refMHChaps/${2}.fasta -o ${3}/RT_${2}
