## Script to assign most appropriate MHC reference haplotype
## to each sample, for assembly scaffolding and ref var calling
## 1/25/23
## wadekj
## Usage: [1] haplotypeRefFile.txt (copied below) [2] genotypes.csv [3] c4genotypes.csv

## MHC Haplotype table
##RefID	OK649233	hg38	OK649231	GL000251.2	OK649235	OK649232	OK649234	OK649236
##MHC Hap	KAS116	PGF	APD	COX	QBL	DBB	MANN	SSTO
##DR	DR1	DR2	DR3	DR3	DR3	DR4	DR4	DR4
##C4 AB	AA	AB	AB	B	A	AB	B	AB
##C4 LS	LL	LL	LL	S	L	LS	L	LL

import sys
from collections import defaultdict
import os
import re
import signal

def main():
    extensionHLA = sys.argv[2].split('.')[1]
    extensionC4 = sys.argv[3].split('.')[1]
    if extensionHLA not in ['csv', 'tsv']: 
        kill_parent_script()
        raise Exception("Invalid HLA genotype format. Please use CSV or TSV file as an input.")
    if extensionC4 not in ['csv', 'tsv']:
        kill_parent_script()
        raise Exception("Invalid C4 genotype format. Please use CSV or TSV file as an input.")
    
    haps=open(sys.argv[1])
    hlaGenos=open(sys.argv[2])
    c4Genos=open(sys.argv[3])
    
    ## hapDict: [DR][C4A/B][C4L/S]=most appropriate ref MHC haplotype
    hapDict=createHapDict(haps)

    ## drGenos: [sampleID][[DRB1A, DRB1B]]
    drGenos=getDRgenos(hlaGenos)

    ## drHaps: [DR allele A, DR allele B]
    drHaps=getDRHaps(drGenos)
    
    ## c4Hap: [sampleID][c4A/B]=c4L/S
    c4Hap=getC4hap(c4Genos)

    assignMHCHaps(hapDict,drHaps,c4Hap)

def createHapDict(haps):
    hapDict=defaultdict(lambda: defaultdict(dict))
    refs=haps.next().strip('\n').split('\t')
    drb=haps.next().strip('\n').split('\t')
    C4AB=haps.next().strip('\n').split('\t')
    C4LS=haps.next().strip('\n').split('\t')
    for i in range(1,len(refs)):
        hapDict[drb[i]][C4AB[i]][C4LS[i]]=refs[i]
    return hapDict

def getDRgenos(hlaGenos):
    drGenos=defaultdict(lambda: defaultdict(list))
    hlaGenos.next()
    fails=[]
    for hla in hlaGenos:
        hla=hla.strip('\n')
        hla=re.split(',|\t', hla)
        sample=hla[0]
        drA=hla[1].split(':')[0]
        drB=hla[2].split(':')[0]
        if drA=='X' or drA=='NA':
            drGenos[sample]['DRB1'].append('X')
        else:
            drGenos[sample]['DRB1'].append(int(drA))
        if drB=='X' or drB=='NA':
            drGenos[sample]['DRB1'].append('X')
        else:
            drGenos[sample]['DRB1'].append(int(drB))

    return drGenos

def getDRHaps(drGenos):
    drHapDict=defaultdict(list)
    for s in sorted(drGenos.keys()):
        drhaps=[]
        for drb in drGenos[s]['DRB1']:
            if drb==1 or drb==10:
                drhaps.append('DR1')
            elif drb==15 or drb==16:
                drhaps.append('DR2')
            elif drb==3 or drb==11 or drb==12 or drb==13 or drb==14:
                drhaps.append('DR3')
            elif drb==4 or drb==9 or drb==7:
                drhaps.append('DR4')
            elif drb==8:
                drhaps.append('DR8')
            elif drb=='X':
                drhaps.append('X')
            else:
                novel=[]
                for g in sorted(drGenos[s].keys()):
                    for u in drGenos[s][g]:
                        novel.append(u)
                novelDR='|'.join(novel)
                drhaps.append(novelDR)
        drHapDict[s]=drhaps
    return drHapDict
            
def getC4hap(c4Genos):
    c4Hap=defaultdict(lambda:defaultdict(str))

    c4Genos.next()
    for c4 in c4Genos:
        c4AB=''
        c4=c4.strip('\n')
        c4=re.split(',|\t', c4)
        sample=c4[7]
        print c4
        c4A=int(c4[1])
        c4B=int(c4[2])
        if c4A >0 and c4B > 0:
            c4AB='AB'
        elif c4A ==1 and c4B== 0:
            c4AB='A'
        elif c4A >1 and c4B == 0:
            c4AB='AA'
        elif c4A==0 and c4B==1:
            c4AB='B'
        elif c4A==0 and c4B>1:
            c4AB='BB'
        elif c4A==0 and c4B==0:
            c4AB='NA'
        else:
            c4AB='CNV'
        c4S=int(c4[4])
        c4L=int(c4[3])
        c4LS=''
        if c4L>0 and c4S>0:
            c4LS='LS'
        elif c4L ==1 and c4S== 0:
            c4LS='L'
        elif c4L >1 and c4S == 0:
            c4LS='LL'
        elif c4L==0 and c4S==1:
            c4LS='S'
        elif c4L==0 and c4S>1:
            c4LS='SS'
        elif c4L==0 and c4S==0:
            c4LS = 'NA'
        else:
            c4LS='CNV'
        c4Hap[sample][c4AB]=c4LS
    return c4Hap

def assignMHCHaps(hapDict,drHap,c4Hap):
    #mmDict=defaultdict(list)
    for sample in sorted(drHap.keys()):
        o=str('./bestHaps/'+sample+'_bestMHChaps.txt')
        out=open(o,'w')
        
        for dr in drHap[sample]:
            try:
                c4AB=c4Hap[sample].keys()[0]
                c4LS=c4Hap[sample][c4AB]
                mhc=hapDict[dr][c4AB][c4LS]
                out.write('%s\t%s\t%s\t%s\t%s\tFull\n'%(sample,mhc,dr, c4AB,c4LS))

            except:
                mhc,mm=handleNovelHaps(hapDict,dr,c4AB,c4LS,sample)
                out.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(sample,mhc,dr, c4AB,c4LS,mm))
    return 

def handleNovelHaps(hapDict,dr,c4AB,c4LS,sample):
    mhc=''
    mm=''
    #novelHap=str(sample+'|'+dr+'|'+c4AB+'|'+c4LS)
    # check for C4ab match
    if hapDict[dr][c4AB]:
        for ls in hapDict[dr][c4AB]:
            if c4LS[0] in ls:
                mhc=hapDict[dr][c4AB][ls]
                #mmDict['DR, C4AB, Partial C4 LS'].append(novelHap)
                mm='DR|C4AB|partC4LS'
            else:
                #mhc=hapDict[dr][c4AB][0]
                for z in hapDict[dr][c4AB]:
                    if hapDict[dr][c4AB] != '':
                        mhc=hapDict[dr][c4AB][z]
                    else:
                        continue
                #mmDict['DR, C4AB'].append(novelHap)
                mm='DR|C4AB'
    # Check for DR mismatch
    elif dr == 'X':
        mhc='X'
        mm='X'
    elif dr not in ['X','DR1','DR2','DR3','DR4','DR8']:
        mhc='Novel'
        #mmDict['Novel'].append(novelHap)
        mm='None'
    # If no C4 match, but DR match, assign first matching DR
    else:
        for d in hapDict[dr].keys():
            for c in hapDict[dr][d].keys():
                if hapDict[dr][d][c]!= '':
                    mhc=hapDict[dr][d][c]
                else:
                    continue
        #mmDict['DR'].append(novelHap)
        mm='DR'
        #mhc=dr
    return mhc,mm

def kill_parent_script():
    try:
        with open("parent_pid.txt", "r") as pid_file:
            parent_pid = int(pid_file.read())
            os.kill(parent_pid, signal.SIGTERM)
    except Exception as e:
        print("Error:", e)

main()
