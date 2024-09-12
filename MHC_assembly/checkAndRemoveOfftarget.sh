#! /bin/bash
## Aligns full supercontigs+contigs from unmapped reads set to full hg38 genome sequence
## Identifies any contigs that exhibit off target mapping elsewhere in the genome and removes them
## Edited: 09/11/24
## Authors: wadekj and susenor


## Args [1] sampleID [2] hg38 [3] assembly file .fa [4] working dir folder </home/LAB_PROJECTS/MHC_PROJECT/4_draftAssemblyAthenav1/sampleID/
## [5] HapB ID [6] HapA ID 

## Map unaligned contigs against hg38 genome
minimap2 -d ../refMHC/${2}.mni ../refMHC/${2}.fa

minimap2 -a ../refMHC/${2}.mni ${3} > ${4}/offTarget/${1}_vs_${2}.sam

minimap2 --cs -c ../refMHC/${2}.fa ${3} > ${4}/offTarget/${1}_UnA_vs_${2}.paf

## If heterozygous, map un-aligned contigs against alternative het haplotype
if [ ${6} != ${5} ]
then 
	minimap2 -d ../refMHC/${5}.mni ../refMHC/${5}.fasta

	minimap2 -a ../refMHC/${5}.mni ${3} > ${4}/offTarget/${1}_vs_${5}_altHap.sam
	minimap2 --cs -c ../refMHC/${5}.fasta ${3} > ${4}/offTarget/${1}_UnA_vs_${5}_altHap.paf

	## Call python script to identify contigs that map outside of MHC
	## Remove contigs with these ids from Amos_supercontigs_unique.fa
	python findAndRemoveOffTarget.py ${4}/offTarget/${1}_UnA_vs_${2}.paf ${4}/offTarget/${1}_UnA_vs_${5}_altHap.paf ${3} ${4}

else
	## Call python script to identify contigs that map outside of MHC
        ## Remove contigs with these ids from Amos_supercontigs_unique.fa
	python findAndRemoveOffTarget.py ${4}/offTarget/${1}_UnA_vs_${2}.paf Null ${3} ${4}
fi

