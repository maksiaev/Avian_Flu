#!/bin/bash
#SBATCH --job-name="H3 2000s Paloma iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s H3_paloma_human_euro_01-01-2000--12-31-2025_aln_trim.fasta -m GTR+G -bb 1000 -nt AUTO


######################
date +"%Y-%m-%d %H:%M"
