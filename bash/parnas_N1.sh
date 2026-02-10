#!/bin/bash
#SBATCH --job-name="NA N1 parnas"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

source myconda
mamba activate base

conda activate parnas_env

parnas -t NA_N1_Paloma_aln_trimmed_deduplicated.fasta.treefile -n 100 --exclude-fully ".*human.*|.*\\|Iberian|.*\\|White" --subtree "NA_N1_Paloma_aln_trimmed_deduplicated_subsampled.treefile"

######################
date +"%Y-%m-%d %H:%M"
