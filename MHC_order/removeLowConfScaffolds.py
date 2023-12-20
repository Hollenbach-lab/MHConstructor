## Script to take in the ragtag.scaffold.confidence.txt output file
## Extract scaffold IDs which have a location confidence <0.3
## Then use those IDs to remove the fasta entry from the assembly file
## 12/5/23

import sys
from collections import defaultdict
from Bio import SeqIO

def main():
    RTdir=sys.argv[1]
    assemblyFasta=sys.argv[2]

    lowQual=QCscaffolds(RTdir)
    removePoorScaff(lowQual,assemblyFasta)

def QCscaffolds(RTdir):
    conf=open(str(RTdir+'/ragtag.scaffold.confidence.txt'))
    lowQual=set()
    conf.next()
    for scaff in conf:
        scaff=scaff.strip('\n').split('\t')
        ## QC2: if float(scaff[2])<=0.1:
        if float(scaff[2])<=0.1:
            print scaff
            lowQual.add(scaff[0])
        elif float(scaff[2])<=0.3 and float(scaff[3])<=0.7:
            print scaff
            lowQual.add(scaff[0])
    return lowQual

def removePoorScaff(lowQual,assembly):
    newOutput=assembly.split('.')[0]
    #newOutput=str(newOutput+'_scaffoldQC2.fasta')
    newOutput=str(newOutput+'_scaffoldQC3.fasta')

    newFasta=open(newOutput,'w')
    with open(assembly) as handle:
        for record in SeqIO.parse(handle,'fasta'):
            fragmentSize=len(record.seq)
            if record.id in lowQual:
                continue
                #print record.id, fragmentSize
            #if record.id in lowQual and fragmentSize< 100000:
                #continue
            else:
                newFasta.write('>%s\n%s\n'% (record.id,record.seq))




main()
