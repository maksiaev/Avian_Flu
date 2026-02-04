#!/bin/bash
#SBATCH --job-name="all_B3.2 alignment"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################

mafft --thread 4 all_B3.2_PB2_2021-11-01--2025-12-05_edit.fasta > all_B3.2_PB2_11-01-2021--12-05-2025_aln.fasta
mafft --thread 4 all_B3.2_PB1_2021-11-01--2025-12-05_edit.fasta > all_B3.2_PB1_11-01-2021--12-05-2025_aln.fasta
mafft --thread 4 all_B3.2_PA_2021-11-01--2025-12-05_edit.fasta > all_B3.2_PA_11-01-2021--12-05-2025_aln.fasta
mafft --thread 4 all_B3.2_HA_2021-11-01--2025-12-05_edit.fasta > all_B3.2_HA_11-01-2021--12-05-2025_aln.fasta
mafft --thread 4 all_B3.2_NP_2021-11-01--2025-12-05_edit.fasta > all_B3.2_NP_11-01-2021--12-05-2025_aln.fasta
mafft --thread 4 all_B3.2_NA_2021-11-01--2025-12-05_edit.fasta > all_B3.2_NA_11-01-2021--12-05-2025_aln.fasta
mafft --thread 4 all_B3.2_MP_2021-11-01--2025-12-05_edit.fasta > all_B3.2_MP_11-01-2021--12-05-2025_aln.fasta
mafft --thread 4 all_B3.2_NS_2021-11-01--2025-12-05_edit.fasta > all_B3.2_NS_11-01-2021--12-05-2025_aln.fasta

######################
date +"%Y-%m-%d %H:%M"
