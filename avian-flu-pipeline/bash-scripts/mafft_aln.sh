#!/bin/bash
set -e
module load mafft

mafft --thread 4 11-01-2021--11-28-2025_concat_573.fasta > 11-01-2021--11-28-2025_concat_573_aln.fasta