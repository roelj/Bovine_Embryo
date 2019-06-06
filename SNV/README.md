# Readme SNPs
SNP information is used to determine the ploidy of cells and to determine the parent-of-origin of genetic events.
These scripts are used to call and filter SNPs in single cell sequencing data and in bulk WGS data.


## Dependencies
1. Opengrid engine
2. bcftools
3. htslib (tabix + bgzip)
4. samtools (1.7)
5. Bos_taurus.UMD3.1.dna.toplevel.fa (reference genome)
6. R (3.5.1)
   - BSgenome.Btaurus.UCSC.bosTau8
   - GenomicRanges

# Readme SNPs
SNP information is used to determine the ploidy of cells and to determine the parent-of-origin of genetic events.
These scripts are used to call and filter SNPs in single cell sequencing data and in bulk WGS data.

## A. Calling and filtering SNPs
#### A1. Merge BAMs by run
- [Project_folder]/Scripts/Other/Merge_BAMs_by_Run.R
- First all single cell bam files are merged by sequencing run into one composite bam file per run. samtools merge *.bam is used to merge all bam files per folder.
```
Rscript [script] [samplesheet_file] [input_folder] [output_folder] [scripts_folder]
eg:
Rscript [Project_folder]/Scripts/Other/Merge_BAMs_by_Run.R [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/BAMs/ [Project_folder]/BAMs/ [Project_folder]/Scripts/
```


#### A2. Merge BAMs from each run into one composite BAM
- [Project_folder]/Scripts/Other/Merge_BAMs_All.sh
- Next the composite bam files of each sequencing run are merged into one large composite bam file containing all reads from all single cells.
Reads previously marked as duplicates are removed and the composite bam file is sorted and indexed. 
```
qsub [Project_folder]/Scripts/Jobs/Merge_BAMs_All.sh
```


#### A3. Call variable SNPs in merged BAM
- [Project_folder]/Scripts/SNV/Call_SNPs_Genome_Merged_Embryo.sh
- This script uses bcftools to call non-ref SNVs in the merged BAM file of the single blastomeres. INDELs are not called (“-I” and “-V indels”). A regions file (“-R”) is used to call only variants on chr1-chr29 and chrX (chrY does not work). ReadGroups are ignored (else variants are called for each read_group (=library) separately). The p-value is set to 0 to call each non-ref variant. The VCF needs post-filtering to remove variants in high coverage regions. 

- Calling command (this is used in the Call_SNPs_Genome_Merged_Embryo.sh script):
```
bcftools mpileup -Ou -f [path_to_reference_genome]/Bos_taurus.UMD3.1.dna.toplevel.fa -R [Project_folder]/Scripts/SNV/Chromosomes_bt8_noY.txt -q10 -Q 20 -d 10000 --ignore-RG [Project_folder]/BAMs/Merged_BAMs/Bovine_Merged_Dedup_Sorted.bam | bcftools call -mO z -P 0 -v -o [Project_folder]/Data/SNVs/SNPs_Merged_Embryo/Bovine_Merged_SNPs_Genome.vcf.gz
```
- How to run:
```
qsub [Project_folder]/Scripts/SNV/Call_SNPs_Genome_Merged_Embryo.sh
```


#### A4. Select heterozygous SNPs
- [Project_folder]/Scripts/SNV/Select_Heterozygous_SNPs.py
- This script selects all heterozygous SNPs in the input VCF. Filters all SNPs from a VCF that have less than [-c N] reads on the alternative or reference allele. It also removes indels and positions that have a read depth of more than [-d N].

```
python [Project_folder]/Scripts/SNV/Select_Heterozygous_SNPs.py -i [input_vcf] -o [output_vcf] -c [minimal_read_counts] -z [zip_output] -d [max_read_depth]

python [Project_folder]/Scripts/SNV/Select_Heterozygous_SNPs.py -i [Project_folder]/Data/SNVs/SNPs_Merged_Embryo/Bovine_Merged_SNPs_Genome.vcf.gz -o /[Project_folder]/Data/SNVs/SNPs_Merged_Embryo/Bovine_Merged_Het_SNPs.vcf -c 2 -z TRUE -d 50
```


#### A5. Genotype sperm DNA per chromosome
- [Project_folder]/Scripts/SNV/Genotype_SNPs_Parallel.R
- The sperm DNA is genotyped for positions of the SNPs that are variable in the composite single cell bam file (from the previous step). The variant calling is parallelized (by chromosome) to speed up the process. Called SNPs per chromosome are merged later in the next step.

```
Rscript [Project_folder]/Scripts/SNV/Genotype_SNPs_Parallel.R [input_VCF] [input_BAM] [input_folder] [output_folder] [output_name] [path_to_reference_genome] 

Example:
Rscript [Project_folder]/Scripts/SNV/Genotype_SNPs_Parallel.R [Project_folder]/Data/SNVs/SNPs_Merged_Embryo/Bovine_Merged_Het_SNPs.vcf [path_to_bovinesperm_bam]/bovinesperm_dedup.realigned.bam [Project_folder]/Scripts/SNV/ [Project_folder]/Data/SNVs/SNPs_Sperm/ BovineSperm_SNPs_Embryo_Het [path_to_reference_genome]/Bos_taurus.UMD3.1.dna.toplevel.fa
```


#### A6. Merge sperm vcfs and select homozygous sperm positions
- [Project_folder]/Scripts/SNV/Filter_BovineSperm_Homozygous_SNPs.sh
- This script first merges the temporary vcfs generated for each chromosome in the previous step into one VCF. This merged VCF contains the genotypes of the sperm for the positions that are variable in the single cells. bcftools is used to filter SNPs (with less than 10 reads or more than 75 reads, and SNPs adjacent to indels) and select only the SNPs homozygous in sperm (output = BovineSperm_Hom_SNPs.vcf.gz). The positions are also exported to bed files.

```
bash [Project_folder]/Scripts/SNV/Filter_BovineSperm_Homozygous_SNPs.sh [Project_folder]/Data/SNVs/SNPs_Sperm/temp/ [Project_folder]/Data/SNVs/SNPs_Sperm/
```


#### A7. Genotype all single cells for homozygous positions in sperm
- [Project_folder]/Scripts/SNV/Call_SNVs_Parallel.R
- bcftools mpileup and bcftools call are used to genotype the single cell bam files for each position homozygous in the father and variable in the composite blastomere bam file.

```
Rscript [Project_folder]/Scripts/SNV/Call_SNVs_Parallel.R TRUE [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Scripts/SNV/ [Project_folder]/BAMs/ [Project_folder]/Data/SNVs/SNPs_Blastomeres_Paternal/ [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_Positions_Sorted.txt.gz [path_to_reference_genome]/Bos_taurus.UMD3.1.dna.toplevel.fa
```


#### A8. Annotate single cells VCFs with paternal genotypes
- [Project_folder]/Scripts/SNV/Annotate_Single_Cell_VCFs.R
- This script is used to annotate each single cell VCF (containing positions that are variable between single cells, but homozygous in the father) with the genotype of the father. 

```
Rscript [Project_folder]/Scripts/SNV/Annotate_Single_Cell_VCFs.R [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Data/SNVs/SNPs_Blastomeres_Paternal/ [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_Positions_Sorted.txt.gz [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_SNPs_header.txt [Project_folder]/Data/SNVs/SNPs_Blastomeres_Paternal_An/ PAT
```

#### A9. Export phased SNP positions
- [Project_folder]/Scripts/SNV/Determine_MAT_SNPs.py
- This scripts determines if the genotype of the single cell is the same or different from the genotype of the father at each position. SNPs different between the father and the blastomere are presumed to be inherited from the mother (labelled with MAT). SNPs overlapping between the father and the blastomere may be paternally or maternally inherited (labelled with PAT/MAT). The phased SNP positions are written to a bed file per single cell.

```
Rscript [Project_folder]/Scripts/SNV/Determine_MAT_SNPs.py -i [Project_folder]/Data/SNVs/SNPs_Blastomeres_Paternal_An/ -o [Project_folder]/Data/SNVs/SNPs_Phased/
```

#### A10. Count MAT and PAT/MAT SNPs per bin
- [Project_folder]/Scripts/SNV/Determine_Maternity.R
- In this step the frequencies of maternal SNPs are calculated per cell and per 10 Mb. This information is later used to determine if cells are uni- or biparental and it is used for visualization. Run both scripts:

```
Rscript [Project_folder]/Scripts/SNV/Determine_Maternity.R cell [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Data/SNVs/SNPs_Phased/ [Project_folder]/Data/SNVs/SNP_Counts/Cell/

Rscript [Project_folder]/Scripts/SNV/Determine_Maternity.R bin [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Data/SNVs/SNPs_Phased/ [Project_folder]/Data/SNVs/SNP_Counts/Bins_10Mb/ bin_size=10000000
```

## B. Strand-specific phasing of single cells
For the cells sequenced by Strand-seq we also want to determine the parent-of-origin of each the DNA strands. The same scripts as in step A are used, but different input is used. First the BAM files of each cell are split into two bam files containing reads of only one of the two strands. Subsequently SNVs are called similarly to normal bam files.

#### B1. Split strand-seq bam files into separate bam files for FW and RV strand
- [Project_folder]/Scripts/SNV/Split_BAMs_Strand.R
- First all single cell bam files are split into separate bam files for the two strands using samtools view -F 20 and samtools view -f 16. These separate bam files are used as input for the other scripts in step B.

```
[Project_folder]/Scripts/SNV/Split_BAMs_Strand.R [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/BAMs/ [Project_folder]/BAMs/
```

#### B2.Genotype both strands for SNP positions homozygous in father
- [Project_folder]/Scripts/SNV/Call_SNVs_Parallel.R
- These scripts can only be run after **step A6** is finished.
```
Rscript [Project_folder]/Scripts/SNV/Call_SNVs_Parallel.R TRUE [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Scripts/SNV/ [Project_folder]/BAMs/BAMs_FW_Strand/ [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_FW/ [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_Positions_Sorted.txt.gz /hpc/cog_bioinf/GENOMES/Bos_taurus_UMC3.1/Bos_taurus.UMD3.1.dna.toplevel.fa Genotype_SNPs_FW

Rscript [Project_folder]/Scripts/SNV/Call_SNVs_Parallel.R TRUE [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Scripts/SNV/ [Project_folder]/BAMs/BAMs_RV_Strand/ [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_RV/ [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_Positions_Sorted.txt.gz /hpc/cog_bioinf/GENOMES/Bos_taurus_UMC3.1/Bos_taurus.UMD3.1.dna.toplevel.fa Genotype_SNPs_RV
```

#### B3. Annotate single cell strand-specific VCFs with paternal genotype
- [Project_folder]/Scripts/SNV/Annotate_Single_Cell_VCFs.R
```
Rscript [Project_folder]/Scripts/SNV/Annotate_Single_Cell_VCFs.R [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_FW/ [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_Positions_Sorted.txt.gz [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_SNPs_header.txt [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_FW_AN/ PAT

Rscript [Project_folder]/Scripts/SNV/Annotate_Single_Cell_VCFs.R [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_RV/ [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_Positions_Sorted.txt.gz [Project_folder]/Data/SNVs/SNPs_Sperm/BovineSperm_Hom_SNPs_header.txt [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_RV_AN/ PAT
```

#### B4. Determine which strand-specific SNPs are maternal or mat/pat
- [Project_folder]/Scripts/SNV/Determine_MAT_SNPs.py
```
python [Project_folder]/Scripts/SNV/Determine_MAT_SNPs.py -i [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_RV_AN/ -o [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_RV_Phased

python [Project_folder]/Scripts/SNV/Determine_MAT_SNPs.py -i [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_FW_AN/ -o [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_FW_Phased
```

#### B5. Count strand-specific MAT and PAT/MAT SNPs per bin
- Rscript [Project_folder]/Scripts/SNV/Determine_Maternity.R
```
Rscript [Project_folder]/Scripts/SNV/Determine_Maternity.R bin [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_FW_Phased/ [Project_folder]/Data/SNVs/SNP_Counts/Bins_10Mb_FW/ bin_size=10000000

Rscript [Project_folder]/Scripts/SNV/Determine_Maternity.R bin [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Data/SNVs/SNPs_Strands/SNPs_RV_Phased/ [Project_folder]/Data/SNVs/SNP_Counts/Bins_10Mb_RV/ bin_size=10000000
```


## C. Genotype single cells to determine ploidy
#### C1. Genotype all single cells for heterozygous positions in embryos
- [Project_folder]/Scripts/SNV/Call_SNVs_Parallel.R
- SNP information is also used to determine the ploidy status of cells. Haploid cells should have no heterozygous variants. Therefore we want to call SNPs in single cells. Just calling each SNP in a single cell gives to many false positives due to technical errors during library prep and sequencing. Therefore we only want to call SNPs are high quality positions. We use genomic positions that were previously determined to be variable in the composite bam (at step A3 and A4). This script can be run after finishing step **A4**. 

```
Rscript [Project_folder]/Scripts/SNV/Call_SNVs_Parallel.R TRUE [Project_folder]/Scripts/Bovine_samplesheet.txt [Project_folder]/Scripts/SNV/ [Project_folder]/BAMs/ [Project_folder]/Data/SNVs/SNPs_Blastomeres/ [Project_folder]/Data/SNVs/SNPs_Sperm/Embryo_Het_SNP_Positions_Sorted.txt.gz [path_to_reference_genome]/Bos_taurus.UMD3.1.dna.toplevel.fa Genotype_blastomeres
```

