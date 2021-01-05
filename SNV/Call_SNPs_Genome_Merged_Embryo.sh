#!/bin/bash

echo "# Start SNV calling"
date

# Path to BAM file
BAM_FILE="/hpc/cuppen/projects/P0017_sc_blastocysts/Ewart_scblastocyst/analysis/roel/SNV/merged_BAMs/total_sorted.bam"

# Output directory
OUTPUT_DIR="/hpc/cuppen/projects/P0017_sc_blastocysts/Ewart_scblastocyst/analysis/roel/SNV/SNV_output"

# Output will be written to this folder:
OUTPUT_FILE="$OUTPUT_DIR/Bovine_Merged_SNPs_Genome.vcf.gz"

# Path to reference
REFERENCE_PATH="/hpc/cog_bioinf/GENOMES/Bos_taurus_UMC3.1/Bos_taurus.UMD3.1.dna.toplevel.fa"

# Create output directory
mkdir -p $OUTPUT_DIR

COMMAND="bcftools mpileup -Ou -f $REFERENCE_PATH -R /hpc/cuppen/projects/P0017_sc_blastocysts/Ewart_scblastocyst/analysis/roel/Bovine_Embryo/SNV/Chromosomes_bt8_noY.txt -q 10 -Q 20 -d 10000 -I --ignore-RG $BAM_FILE | bcftools call -mO z -P 0 -v -V indels -o $OUTPUT_FILE"
echo "$COMMAND"
eval "$COMMAND"

COMMAND="bcftools index $OUTPUT_FILE"
echo "$COMMAND"
eval "$COMMAND"

COMMAND="bcftools stats $OUTPUT_FILE > $OUTPUT_FILE.vchk"
echo "$COMMAND"
eval "$COMMAND"

echo "# Finished SNV calling"
date
