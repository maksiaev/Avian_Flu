#!/bin/bash
#SBATCH --job-name="A3 iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

module load iqtree
srun  iqtree2 -s A3_full_genome_11-01-2021--06-13-2025_aln.fas -m GTR+G -bb 1000 -nt AUTO

######################
date +"%Y-%m-%d %H:%M"
