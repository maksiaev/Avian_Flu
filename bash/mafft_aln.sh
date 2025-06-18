#!/bin/bash
set -e
module load mafft

mafft --thread 4 B3_13_PB1_concat_11-01-2021--06-13-2025.fasta > B3_13_PB1_concat_11-01-2021--06-13-2025_aln.fasta