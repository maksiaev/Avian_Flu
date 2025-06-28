#!/bin/bash
#SBATCH --job-name="all D1.1 iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alvin.crespobellido@nih.gov"
##################################################
######################

module load iqtree
srun  iqtree2 -s newid-D1.1_concat_genome_MAY14_trim_aln_n2240_FINAL.fasta -m GTR+G -bb 1000 -nt AUTO

######################
date +"%Y-%m-%d %H:%M"
