#! /bin/bash
## Script to blast full MHC assembly to a ref seq
## Extract variation from the output blast files
## convert into vcf format
## 1/20/23

## updated: 08/18/23
## Aligns full supercontig set to the ref, instead of just unassembled
## wadekj


## Args [1] sampleID [2] hg38 [3] assembly file .fa [4] working dir folder </home/LAB_PROJECTS/MHC_PROJECT/4_draftAssemblyAthenav1/sampleID/
## [5] HapB ID [6] HapA ID 
## Map unaligned contigs against hg38 genome
minimap2 -d ../refMHChaps/${2}.mni ../refMHChaps/${2}.fa

minimap2 -a ../refMHChaps/${2}.mni ${3} > ${4}/offTarget/${1}_vs_${2}.sam

minimap2 --cs -c ../refMHChaps/${2}.fa ${3} > ${4}/offTarget/${1}_UnA_vs_${2}.paf

## If heterozygous, map un-aligned contigs against alternative het haplotype
if [ ${6} != ${5} ]
then 
	minimap2 -d ../refMHChaps/${5}.mni ../refMHChaps/${5}.fasta

	minimap2 -a ../refMHChaps/${5}.mni ${3} > ${4}/offTarget/${1}_vs_${5}_altHap.sam
	minimap2 --cs -c ../refMHChaps/${5}.fasta ${3} > ${4}/offTarget/${1}_UnA_vs_${5}_altHap.paf

	## Call python script to identify contigs that map outside of MHC
	## Remove contigs with these ids from Amos_supercontigs_unique.fa
	python findAndRemoveOffTarget.py ${4}/offTarget/${1}_UnA_vs_${2}.paf ${4}/offTarget/${1}_UnA_vs_${5}_altHap.paf ${3} ${4}

else
	## Call python script to identify contigs that map outside of MHC
        ## Remove contigs with these ids from Amos_supercontigs_unique.fa
	python findAndRemoveOffTarget.py ${4}/offTarget/${1}_UnA_vs_${2}.paf Null ${3} ${4}
fi

