#!/bin/bash
#SBATCH --job-name="MAFFT cat genome"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="alexander.maksiaev@nih.gov"
##################################################
######################


set -e
module load mafft

mafft --thread $SLURM_CPUS_PER_TASK --reorder --add B3.13_genome_combined_06-10-2025_clean.fasta newid-B3.13_concat_genome_MAY14_trim_aln_n3651_FINAL_updated_06-11-2025_clean.fasta > newid-B3.13_genome_JUNE6_aln.fasta 

for file_name in $(ls); do mafft --thread $SLURM_CPUS_PER_TASK --reorder --add B3.13_genome_combined_06-10-2025_clean.fasta newid-B3.13_concat_genome_MAY14_trim_aln_n3651_FINAL_updated_06-11-2025_clean.fasta done > newid-B3.13_genome_JUNE6_aln.fasta 