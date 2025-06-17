#!/bin/bash
#SBATCH --job-name="B3_13 MAFFT concatenation genome"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

set -e
module load mafft

segments[0]="PB2"
segments[1]="PB1"
segments[2]="PA"
segments[3]="HA"
segments[4]="NP"
segments[5]="NA"
segments[6]="MP"
segments[7]="NS"

for segment in ${segments[@]}; do mafft --thread --cpus-per-task=4 B3_13_${segment}_concat_11-01-2021--06-13-2025.fasta > B3_13_${segment}_concat_11-01-2021--06-13-2025_aln.fasta; done

######################
date +"%Y-%m-%d %H:%M"