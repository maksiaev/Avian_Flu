#!/bin/bash
#SBATCH --job-name="NS Paloma iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s NS_Paloma_aln_trimmed_deduplicated.fasta -m GTR+G -bb 1000 -nt AUTO

######################
date +"%Y-%m-%d %H:%M"
