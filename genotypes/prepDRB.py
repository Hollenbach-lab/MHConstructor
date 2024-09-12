import os
import pandas as pd
import sys

# A supporting script for genotypeDRB.sh to format the HLA-DRB1 genotypes to a nice tabular format
## Edited: 9/11/24
## Authors: susenor and wadekj

def main():
    # alleles = "/home/rsuseno/MHConstructor/genotypes/alleles"
    alleles = sys.argv[1]
    directory = os.fsencode(alleles)
    drb11 = []
    drb12 = []
    sid = []

    for file in os.listdir(directory):
        file = file.decode('utf-8')
        sid.append(file.split('_')[0])
        toread = f'{alleles}/{file}'
        allele = pd.read_csv(toread, sep = '\t', header = None)
        tmp1=(allele.iloc[0,0].split(' ')[0].split('*')[1])

        # Check whether there's only a single genotype in the _allele.tsv file
        if len(allele) == 1:
            tmp2=tmp1
        else:
            tmp2=(allele.iloc[1,0].split(' ')[0].split('*')[1])

        tmp1 = tmp1.split(':')[0]
        tmp2 = tmp2.split(':')[0]

        # Uncomment for second field resolution
        # Limit the result to only include 2 resolutions: cutting it off if resolution>2
        # test1 = len(tmp1.split(':'))
        # test2 = len(tmp2.split(':'))
        # if test1 > 2:
        #     spl1 = tmp1.split(':')
        #     tmp1 = spl1[0]+':'+spl1[1]
        # if test2 > 2:
        #     spl2 = tmp2.split(':')
        #     tmp2 = spl2[0]+':'+spl2[1]
        
        drb11.append(tmp1)
        drb12.append(tmp2)

    final_df = pd.DataFrame()
    final_df['SampleID'] = sid
    final_df['HLA-DRB1.1'] = drb11
    final_df['HLA-DRB1.2'] = drb12
    final_df.to_csv('HLAgenotypes.csv', index=False)

main()
