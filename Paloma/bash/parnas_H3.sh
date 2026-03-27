#!/bin/bash
#SBATCH --job-name="HA H3 parnas"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

source myconda
mamba activate base

conda activate parnas_env

parnas -t HA_H3_Paloma_aln_trimmed_deduplicated.fasta.treefile -n 750 --exclude-fully ".*human.*|.*\\|Iberian|.*\\|White" --subtree "HA_H3_Paloma_aln_trimmed_deduplicated_subsampled.treefile"

######################
date +"%Y-%m-%d %H:%M"
