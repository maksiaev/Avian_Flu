#!/bin/bash
#SBATCH --job-name="HA H3 Paloma iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s N2_paloma_human_euro_h1n1_01-01-1976--12-31-2025_aln_trim.fasta -m GTR+G -bb 1000 -nt AUTO


######################
date +"%Y-%m-%d %H:%M"
