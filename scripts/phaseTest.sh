#!/bin/bash

# Region to run script on
region=$1
options=$2

# Record commands
set -o xtrace

# Run margin
echo Running Margin
time ./margin ../tests/shasta_phasing_data.100kb_5x/${region}/*.bam ../tests/shasta_phasing_data.100kb_5x/${region}/HG002.shasta.*.fasta ../params/allParams.np.json ${options}  --logLevel DEBUG

# Calculate identity
echo Comparing to haplotype1
time python3 ../scripts/dirty_assembly_compare.py ../tests/shasta_phasing_data.100kb_5x/${region}/HG002_h1*.fa output.fa

echo Comparing to haplotype2
time python3 ../scripts/dirty_assembly_compare.py ../tests/shasta_phasing_data.100kb_5x/${region}/HG002_h2*.fa output.fa

# Fasta split the phased output
python3 ../scripts/fasta_split.py output.fa

# Build collection of differences
time python3 ../scripts/dirty_assembly_compare.py output.fa_0 output.fa_1 verbose > predictedMismatches.txt
time python3 ../scripts/dirty_assembly_compare.py ../tests/shasta_phasing_data.100kb_5x/${region}/HG002_h1*.fa ../tests/shasta_phasing_data.100kb_5x/${region}/HG002_h2*.fa verbose > trueMismatches.txt

# Compare differences
echo Comparing predicted hets
time python3 ../scripts/compareHets.py trueMismatches.txt predictedMismatches.txt


