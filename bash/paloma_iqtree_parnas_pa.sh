#!/bin/bash
#SBATCH --job-name="Paloma iqtree PA post parnas"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s PA_parnas_humans_Paloma_aln_trimmed.fasta -m GTR+G -bb 1000 -nt AUTO

######################
date +"%Y-%m-%d %H:%M"
