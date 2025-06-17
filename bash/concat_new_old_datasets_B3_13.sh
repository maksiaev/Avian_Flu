#!/bin/bash
#SBATCH --job-name="B3.13 old and new concatenation"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

segments[0]="PB2"
segments[1]="PB1"
segments[2]="PA"
segments[3]="HA"
segments[4]="NP"
segments[5]="NA"
segments[6]="MP"
segments[7]="NS"

for segment in ${segments[@]}; do cat *"$segment"*.fasta > B3_13_${segment}_concat_11-01-2021--06-13-2025.fasta; done

######################
date +"%Y-%m-%d %H:%M"