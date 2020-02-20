#!/bin/bash

# Region to run script on
region=$1
depth=$2
rle=$3

# Record commands
set -o xtrace

# Run margin
echo Running Margin
time ./margin ../tests/shasta_phasing_data.100kb_5x/${region}/*${depth}x.bam ../tests/shasta_phasing_data.100kb_5x/${region}/*shasta.fasta ../params/allParams.np${rle}.json  --logLevel DEBUG

# Calculate identity
echo Comparing to truth
time python3 ../scripts/dirty_assembly_compare.py ../tests/shasta_phasing_data.100kb_5x/${region}/*truth.fasta output.fa
