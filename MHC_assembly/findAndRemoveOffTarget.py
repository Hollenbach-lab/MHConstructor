## Script to process .paf output and find unassembled contigs that map
## outside MHC
## 03/20/23
## wadekj

import sys
from collections import defaultdict
from Bio import SeqIO
def main():
    paf=open(sys.argv[1])
    altHetPaf=sys.argv[2]
    unAcontigs=sys.argv[3]
    unAdir=sys.argv[4]

    ## filter out hg38 genome-mapping contigs
    offTarget=filterOffTarget(paf)
    removeOffTarget(offTarget,unAcontigs,unAdir)

    ## if heterozygous individual, filter out HapB-mapping contigs 
    if altHetPaf != 'Null':
        altHetPaf=open(sys.argv[2])
        offTargetHet=filterHetOT(altHetPaf)
        offTarget=offTarget.union(offTargetHet)

    ## Write new fasta without off-target mapping contigs
    removeOffTarget(offTarget,unAcontigs,unAdir)

    
def filterOffTarget(paf):
    offTarget=set()
    for line in paf:
        data=line.strip('\n').split('\t')
        contig=data[0]
        mapChrom=data[5]
        mapStart=int(data[7])
        mapEnd=int(data[8])
        if 'chr6' in mapChrom:
            if 'alt' in mapChrom:
                continue
            elif mapStart >= 28525013 or mapStart<= 33457522:
                continue
            elif mapEnd >= 28525013 or mapEnd <= 33457522:
                continue
            else:
                offTarget.add(contig)
        else:
            offTarget.add(contig)
    return offTarget

def filterHetOT(paf):
    offTarget=set()
    for line in paf:
        data=line.strip('\n').split('\t')
        contig=data[0]
        offTarget.add(contig)
    return offTarget
    
def removeOffTarget(offTarget,contigs,unAdir):
    out=open(str(unAdir+'/Unassembled_velvet_OT.fa'),'w')
    for record in SeqIO.parse(contigs, "fasta"):
        if record.id in offTarget:
            continue
        else:
            out.write('>%s\n%s\n'%(record.id, record.seq))

main()
