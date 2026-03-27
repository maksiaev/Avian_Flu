#!/bin/bash
#SBATCH --job-name="HA H3 Paloma iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s subsampled_europe_human+swine_pdmOnly_2009-2026_aln_trim_downsampled.fasta -m GTR+G -bb 1000 -nt AUTO


######################
date +"%Y-%m-%d %H:%M"
