## script to rename imgt alleles so they are HLA alleles, not ids
## 01/03/24
## wadekj

import sys
from Bio import SeqIO

def main():
    f=open(sys.argv[1])
    out=open('DRB1_IMGT_010324.fasta','w')
    for record in f:
        if record.startswith('>'):
            name=record.split(' ')[1]
            out.write('>%s\n'%(name))
        else:
            out.write('%s\n'%(record.strip('\n')))

main()

