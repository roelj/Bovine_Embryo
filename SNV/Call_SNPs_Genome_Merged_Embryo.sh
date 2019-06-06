#!/bin/bash

#$ -S /bin/bash
#$ -l h_vmem=30G
#$ -l h_rt=08:00:00
#$ -cwd
#$ -o /hpc/cog_bioinf/cuppen/project_data/Ewart_Single_Cell/Bovine/Scripts/SNV/Logs/20190527_Call_SNPs_Genome_Merged_Embryo_log.txt
#$ -e /hpc/cog_bioinf/cuppen/project_data/Ewart_Single_Cell/Bovine/Scripts/SNV/Logs/20190527_Call_SNPs_Genome_Merged_Embryo_log.txt

echo "# Start SNV calling"
date

guixr load-profile ~/.guix-profile/ -- <<EOF
	# Path to BAM file
	BAM_FILE=/hpc/cog_bioinf/cuppen/project_data/Ewart_Single_Cell/Bovine/BAMs/Merged_BAMs/Bovine_Merged_Dedup_Sorted.bam

	# Output will be written to this folder:
	OUTPUT_FILE=/hpc/cog_bioinf/cuppen/project_data/Ewart_Single_Cell/Bovine/Data/SNVs/SNPs_Merged_Embryo/Bovine_Merged_SNPs_Genome.vcf.gz
	
	# Path to reference
	REFERENCE_PATH=/hpc/cog_bioinf/GENOMES/Bos_taurus_UMC3.1/Bos_taurus.UMD3.1.dna.toplevel.fa

	# Create output directory
	mkdir -p /hpc/cog_bioinf/cuppen/project_data/Ewart_Single_Cell/Bovine/Data/SNVs/SNPs_Merged_Embryo/

	echo "bcftools mpileup -Ou -f \$REFERENCE_PATH -R /hpc/cog_bioinf/cuppen/project_data/Ewart_Single_Cell/Bovine/Scripts/SNV/Chromosomes_bt8_noY.txt -q 10 -Q 20 -d 10000 --ignore-RG \$BAM_FILE | bcftools call -mO z -P 0 -v -o \$OUTPUT_FILE"
	bcftools mpileup -Ou -f \$REFERENCE_PATH -R /hpc/cog_bioinf/cuppen/project_data/Ewart_Single_Cell/Bovine/Scripts/SNV/Chromosomes_bt8_noY.txt -q 10 -Q 20 -d 10000 -I --ignore-RG \$BAM_FILE | bcftools call -mO z -P 0 -v -V indels -o \$OUTPUT_FILE
	
	echo "bcftools index \$OUTPUT_FILE"
	bcftools index \$OUTPUT_FILE
	
	echo "bcftools stats \$OUTPUT_FILE > \$OUTPUT_FILE.vchk"
	bcftools stats \$OUTPUT_FILE > \$OUTPUT_FILE.vchk	

EOF

echo "# Finished SNV calling"
date
