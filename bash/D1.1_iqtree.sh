#!/bin/bash
#SBATCH --job-name="all D1.1 iqtree"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alvin.crespobellido@nih.gov"
##################################################
######################

module load iqtree
iqtree2 -s D1_1_11-01-2021--07-18-2025_concat_2905.fasta -m GTR+G -bb 1000 -nt AUTO

######################
date +"%Y-%m-%d %H:%M"
