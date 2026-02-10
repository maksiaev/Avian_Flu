#!/bin/bash
#SBATCH --job-name="MP parnas"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

source myconda
mamba activate base

conda activate parnas_env

parnas -t MP_Paloma_aln_trimmed_deduplicated.fasta.treefile -n 100 --exclude-fully ".*human.*|.*\\|Iberian|.*\\|White" --subtree "MP_Paloma_aln_trimmed_deduplicated_subsampled.treefile"

######################
date +"%Y-%m-%d %H:%M"
