#!/bin/bash
#SBATCH --job-name="MAFFT cat genome"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################


set -e
module load mafft

mafft --thread $SLURM_CPUS_PER_TASK --reorder --add B3_13_concat_07-05-2025--07-18-2025.fas B3_13_11-01-2021--07-04-2025_concat_4819_updated_07-21-2025.fasta > B3_13_11-01-2021--07-18-2025_concat.fasta 

for file_name in $(ls); do mafft --thread $SLURM_CPUS_PER_TASK --reorder --add D1_1_concat_07-05-2025--07-18-2025.fas D1_1_11-01-2021--07-04-2025_concat_2620_updated_07-21-2025.fasta done > D1_1_11-01-2021--07-18-2025_concat.fasta 