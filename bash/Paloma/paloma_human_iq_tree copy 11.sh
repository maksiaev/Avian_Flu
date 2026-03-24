#!/bin/bash
#SBATCH --job-name="HA H3 Paloma iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s PB1_Paloma_wLineage_wDatesAdded_v7_aln_trim.fasta -m GTR+G -bb 1000 -nt AUTO


######################
date +"%Y-%m-%d %H:%M"
