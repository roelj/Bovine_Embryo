#!/bin/bash

# run in qlogin
# Argument 1 = temp_folder containing the SNPs per chromosome
# Argument 2 = output_folder

# This script first merges the temporary vcfs generated for each chromosome by "Genotype_SNPs_Parallel.R" into one VCF. This VCF contains the sperm SNP info about the positions heterozygous in the embryos.
echo "# Merging VCFs"

COMMAND="bcftools concat $1*.vcf.gz -O z -o $2BovineSperm_SNPs_Het_Embryo.vcf.gz"
echo "$COMMAND"
eval "$COMMAND"

COMMAND="bcftools sort -T $(pwd) $2BovineSperm_SNPs_Het_Embryo.vcf.gz -O z"
echo "$COMMAND"
eval "$COMMAND"

COMMAND="bcftools index $2BovineSperm_SNPs_Het_Embryo.vcf.gz"
echo "$COMMAND"
eval "$COMMAND"

COMMAND="bcftools stats $2BovineSperm_SNPs_Het_Embryo.vcf.gz > $2BovineSperm_SNPs_Het_Embryo.stats"
echo "$COMMAND"
eval "$COMMAND"

echo "# Select high quality SNV calls"

# First remove SNPs next to indels
bcftools norm -d all $2BovineSperm_SNPs_Het_Embryo.vcf.gz -O z -o $2Deduplicated_positions.vcf.gz
bcftools index $2Deduplicated_positions.vcf.gz

# select the positions that have to be removed
bcftools isec -C $2BovineSperm_SNPs_Het_Embryo.vcf.gz $2Deduplicated_positions.vcf.gz -O z -o $2Duplicated_positions.vcf.gz -w1

bcftools view -T ^$2Duplicated_positions.vcf.gz $2BovineSperm_SNPs_Het_Embryo.vcf.gz -O z -o $2BovineSperm_SNPs_Het_Embryo_Filtered.vcf.gz
bcftools index $2BovineSperm_SNPs_Het_Embryo_Filtered.vcf.gz

# Remove variant calls with a depth lower than 10 or higher than 75. Also remove indels
bcftools view -e "DP>75 | DP<10 | IDV>0" -V indels $2BovineSperm_SNPs_Het_Embryo_Filtered.vcf.gz -o $2BovineSperm_SNPs_Het_Embryo_Filtered2.vcf
bgzip -c $2BovineSperm_SNPs_Het_Embryo_Filtered2.vcf > $2BovineSperm_SNPs_Het_Embryo_Filtered2.vcf.gz
bcftools index $2BovineSperm_SNPs_Het_Embryo_Filtered2.vcf.gz
bcftools stats $2BovineSperm_SNPs_Het_Embryo_Filtered2.vcf.gz > $2BovineSperm_SNPs_Het_Embryo_Filtered2.stats

echo "# Selecting homozygous SNPs"
bcftools view -g hom -V indels $2BovineSperm_SNPs_Het_Embryo_Filtered2.vcf.gz > $2BovineSperm_Hom_SNPs.vcf

echo "# Zip and index vcfs"
bgzip -c $2BovineSperm_Hom_SNPs.vcf > $2BovineSperm_Hom_SNPs.vcf.gz
bcftools index $2BovineSperm_Hom_SNPs.vcf.gz
bcftools stats $2BovineSperm_Hom_SNPs.vcf.gz > $2BovineSperm_Hom_SNPs.stats
echo "# Zip and index vcfs"


# Write the homozygous position to a bed file (chr start end GT (=10th column)
echo "# Select homozygous SNP positions with paternal genotype"
echo "sed -e 's/chr//' $2BovineSperm_Hom_SNPs.vcf | awk '{OFS=\"\t\"; if (!/^#/){print $1,$2,$2,$10}}' > $2BovineSperm_Hom_Positions.txt"

sed -e 's/chr//' $2"BovineSperm_Hom_SNPs.vcf" | awk '{OFS="\t"; if (!/^#/){print $1,$2,$2,substr($10,1,3)}}' > $2"BovineSperm_Hom_Positions.txt"

echo "# Sort, zip and index SNP positions file. Remove duplicated lines."
sort -k1,1 -k2,2n $2"BovineSperm_Hom_Positions.txt" | uniq -u > $2"BovineSperm_Hom_Positions_Sorted.txt"
# Set header (necessary for later annotations)
sed -i $'1 i\\\n#CHROM\tFROM\tTO\tPAT' $2"BovineSperm_Hom_Positions_Sorted.txt"

# Generate file containing header for annotating VCF
echo "##INFO=<ID=PAT,Number=1,Type=Integer,Description=\"Genotype of father\">" > $2"BovineSperm_Hom_SNPs_header.txt"


# Write the heterozygous positions in the embryo to a bed file (chr start end GT (=10th column)
echo "# Select all SNP positions from the embryonic composite VCF"
echo "sed -e 's/chr//' $2BovineSperm_SNPs_Het_Embryo_Filtered2.vcf | awk '{OFS=\"\t\"; if (!/^#/){print $1,$2,$2}}' > $2Embryo_Het_SNP_Positions.txt"

sed -e 's/chr//' $2"BovineSperm_SNPs_Het_Embryo_Filtered2.vcf" | awk '{OFS="\t"; if (!/^#/){print $1,$2,$2}}' > $2"Embryo_Het_SNP_Positions.txt"

echo "# Sort, zip and index SNP positions file. Remove duplicated lines."
sort -k1,1 -k2,2n $2"Embryo_Het_SNP_Positions.txt" | uniq -u > $2"Embryo_Het_SNP_Positions_Sorted.txt"

# zip and index the SNP position files
bgzip $2"BovineSperm_Hom_Positions_Sorted.txt"
tabix -s1 -b2 -e3 $2"BovineSperm_Hom_Positions_Sorted.txt.gz"

bgzip -c $2"Embryo_Het_SNP_Positions_Sorted.txt" > $2"Embryo_Het_SNP_Positions_Sorted.txt.gz"
tabix -s1 -b2 -e3 $2"Embryo_Het_SNP_Positions_Sorted.txt.gz"

echo "Finished"


