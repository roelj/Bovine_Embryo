# This script is used to determine if cells are haploid, uniparental, diploid or multiploid based on copy number states and SNP distribution

options(scipen = 999)
args <- commandArgs(trailingOnly=TRUE)
samplesheet_file <- args[1] # [project_folder]/Data/QC/Merged_quality_metrics.txt
input_folder <- args[2] # should be [project_folder]/Data/
results_folder <- args[3] # plots will be written in this folder

print("# Determining ploidy of single cells")
library(ggplot2)

print ("# Reading data")
## Read all input data
if(file.exists(samplesheet_file) == TRUE){
  samplesheet <- read.delim(samplesheet_file, stringsAsFactors = F)
} else {
  stop("! Samplesheet not found")
}

resolution <- 1000000
QC_folder <- paste(input_folder, "QC/", sep = "")
StrandSeq_metrics <- read.delim(paste(QC_folder, "/Merged_StrandSeq_metrics.txt", sep = ""), stringsAsFactors = F)
Merged_idxstats <- read.delim(paste(QC_folder, "/Merged_idxstats.txt", sep = ""), stringsAsFactors = F, check.names = F)
SNP_Stats <-  read.delim(paste(QC_folder, "Merged_SNP_Stats.txt",sep = ""), stringsAsFactors = F)
MAT_SNPs <- read.delim(paste(input_folder, "SNVs/SNP_Counts/Cell/PAT_MAT_SNPs_per_cell.txt", sep = ""), stringsAsFactors = F)

if(dir.exists(results_folder) == FALSE){
  dir.create(results_folder, recursive = T)
}

## SNPS
## Strand-seq
## YX (absence of X)
## Deletions are nullosomies 

## Determine the number of sex chromosomes in the cell based on reads on X/Y chromosomes
XY_Overview <- samplesheet[,c("Cell", "Embryo","Stage","Group","Method", "total.read.count")]

# The total.read.count is determined by Aneufinder, but it excludes reads overlapping blacklisted regions
Reads_per_cell <- data.frame(Cell = names(colSums(Merged_idxstats[,-c(1,2)])), Reads = colSums(Merged_idxstats[,-c(1,2)]))

y_reads <- data.frame(Cell = names(Merged_idxstats)[-1], 
                      X_reads = as.numeric(Merged_idxstats[Merged_idxstats[,1] == "X", -1]),
                      Y_reads = as.numeric(Merged_idxstats[Merged_idxstats[,1] == "Y", -1]), stringsAsFactors = F)
y_reads <- y_reads[which(y_reads$Cell != "V2"),] # Remove chromosome lengths
XY_Overview <- merge(XY_Overview, y_reads, by = "Cell")
XY_Overview <- merge(XY_Overview, Reads_per_cell, by = "Cell")

# Portion of reads mapping to X or Y chromosomes:
XY_Overview$Y_reads_rel <- XY_Overview$Y_reads / XY_Overview$Reads
XY_Overview$X_reads_rel <- XY_Overview$X_reads / XY_Overview$Reads
XY_Overview$YX <- XY_Overview$Y_reads/XY_Overview$X_reads

XY_Overview$Sex <- ifelse(XY_Overview$Y_reads_rel > 0.0015 & XY_Overview$X_reads_rel > 0.015 &
                            XY_Overview$Reads > 50000, "XY", "ND")

XY_Overview$Sex <- ifelse(XY_Overview$Y_reads_rel > 0.0015 & XY_Overview$X_reads_rel < 0.015 &
                            XY_Overview$Reads > 50000, "Y", XY_Overview$Sex )

XY_Overview$Sex <- ifelse(XY_Overview$Y_reads_rel < 0.0015 & XY_Overview$X_reads_rel > 0.025 &
                            XY_Overview$Reads > 50000, "XX", XY_Overview$Sex )

summary(factor(XY_Overview$Sex))

# subsequently determine the sexes of the embryos by looking at majority
for(embryo in unique(XY_Overview$Embryo)){
  
  embryo_data <- XY_Overview[XY_Overview$Embryo == embryo,]
  
  male_cells <- length(which(embryo_data$Sex == "XY" | embryo_data$Sex == "Y"))
  female_cells <- length(which(embryo_data$Sex == "XX"))
  
  if((male_cells+female_cells) / nrow(embryo_data) > 0.5){
    XY_Overview$Embryo_Sex[XY_Overview$Embryo == embryo] <- ifelse(male_cells/(female_cells+0.1) > 1.5, "XY", "ND")
    XY_Overview$Embryo_Sex[XY_Overview$Embryo == embryo] <- ifelse(female_cells/(male_cells+0.1) > 1.5, "XX", 
                                                                      XY_Overview$Embryo_Sex[XY_Overview$Embryo == embryo])
  } else {
    XY_Overview$Embryo_Sex[XY_Overview$Embryo == embryo] <- "ND"
  }
}

write.table(x = XY_Overview, file = paste(QC_folder, "Overview_chrXY.txt", sep = ""), quote = F, sep = "\t")

# Plot overviews of reads on X/Y chromosomes
pdf(file = paste(results_folder, "XY_analysis.pdf", sep = ""), width = 4, height = 3, pointsize = 10)
ggplot(XY_Overview[XY_Overview$total.read.count > 100000,], aes(x = Y_reads_rel*100, y = X_reads_rel*100, col = Sex)) + geom_point() +
  geom_vline(xintercept = 0.0015*100) +
  geom_hline(yintercept = c(0.015*100)) + labs(x = "Reads on chrY (% of total)", y = "Reads on xhrX (% of total)") + scale_color_brewer(palette="Set1") + theme_classic()
ggplot(XY_Overview[XY_Overview$total.read.count > 100000,], aes(x = Y_reads_rel*100, fill = Sex)) + geom_density(alpha = 0.7)+ 
    labs(x = "Reads on chrY (% of total)")+ scale_fill_brewer(palette="Set1")+ theme_classic()
ggplot(XY_Overview[XY_Overview$total.read.count > 100000,], aes(x = X_reads_rel*100, fill = Sex)) + geom_density(alpha = 0.7) + 
  labs(x = "Reads on chrX (% of total)")+ scale_fill_brewer(palette="Set1")+ theme_classic()
ggplot(XY_Overview[XY_Overview$total.read.count > 100000,], aes(x = Stage, fill =  Sex)) + geom_bar(stat = "count", position = "dodge", col = "black") + 
  scale_fill_brewer(palette="Set1")+ theme_classic() + scale_y_continuous(expand = c(0,0)) + labs(y = "Cells (#)")

ggplot(XY_Overview[XY_Overview$total.read.count > 100000,], aes(x = Stage, fill = Embryo_Sex)) + geom_bar(stat = "count", position = "dodge", col = "black") + 
  scale_fill_brewer(palette="Set1")+ theme_classic() + scale_y_continuous(expand = c(0,0)) + labs(y = "Embryos (#)")

dev.off()

## 
Ploidy_overview <- samplesheet[,c("Cell", "Embryo","Stage","Group","Method", "Run", "total.read.count")]

## Ploidy based on heterozygous SNPs
SNP_Stats$AF0AF1 <- SNP_Stats$AF0_Filtered / (SNP_Stats$AF1/1000)
SNP_Stats$Ploidy_SNPs <- ifelse(SNP_Stats$AF0AF1 < 1 & SNP_Stats$Total_SNPs > 3000, "Haploid/Uniparental", "ND")
SNP_Stats$Ploidy_SNPs <- ifelse(SNP_Stats$AF0AF1 > 1 & SNP_Stats$Total_SNPs > 3000, "Diploid", SNP_Stats$Ploidy_SNPs)

Ploidy_overview <- merge(Ploidy_overview, SNP_Stats[,c("Cell", "Total_SNPs", "AF0AF1", "AF0", "AF0_Filtered","Ploidy_SNPs")])

## Ploidy based on absence of X-chr
Ploidy_XY <- XY_Overview[,c("Cell", "Y_reads_rel", "X_reads_rel", "YX", "Sex")]
Ploidy_XY$Ploidy_XY <- ifelse(Ploidy_XY$Sex == "Y", "Haploid", "ND")
Ploidy_XY$Ploidy_XY <- ifelse(Ploidy_XY$Sex == "XY", "Diploid", Ploidy_XY$Ploidy_XY)

Ploidy_overview <- merge(Ploidy_overview, Ploidy_XY, by = "Cell", all.x = T)

## Deletions in haploid cells should be nullisomies and therefore should not contain reads. 
# Check the amount of reads per deletion in each cell
Deletion_Overview <- data.frame()

CNV_folder <- paste(input_folder,"CNVs_cells/", resolution, "/", sep  = "")
for(Run in unique(Ploidy_overview$Run[Ploidy_overview$Method == "WGS" | Ploidy_overview$Stage == "2-cell"])){
  print(Run)
  for(Cell in Ploidy_overview$Cell[Ploidy_overview$Run == Run]){
    CNVs_File_Cell <- paste(CNV_folder, Run, "/CNVs_", Cell, ".txt", sep = "")
    if(file.exists(CNVs_File_Cell)){
      CNVs_Cell <- read.delim(CNVs_File_Cell)
      Deletions_Cell <- CNVs_Cell[which(CNVs_Cell[,Cell] < 0 & CNVs_Cell$width > 10e6 & CNVs_Cell$seqnames != "X"),]
      if(nrow(Deletions_Cell) > 0){
        Deletion_Overview_Cell <- data.frame(Cell = Cell, Deletions = nrow(Deletions_Cell), 
                                             Nullisomies = nrow(Deletions_Cell[which(Deletions_Cell$Relative_CNV_Change < 0.1),]))
        Deletion_Overview <- rbind(Deletion_Overview, Deletion_Overview_Cell)
      }
    }
  }
}

Deletion_Overview$Null_Dels <- Deletion_Overview$Nullisomies / Deletion_Overview$Deletions
Deletion_Overview$DEL0 <- ifelse(Deletion_Overview$Null_Dels > 0.5, "Partial", "No")
Deletion_Overview$DEL0 <- ifelse(Deletion_Overview$Null_Dels > 0.8, "Yes", Deletion_Overview$DEL0)

Ploidy_overview <- merge(Ploidy_overview, Deletion_Overview[,c("Cell", "DEL0")], all.x = T)
Ploidy_overview$DEL0[is.na(Ploidy_overview$DEL0)] <- "ND"

# Haploid and uniparental cells should have a bias in the SNPs inherited from the mother or the father 

Ploidy_overview <- merge(Ploidy_overview, MAT_SNPs[,c("Cell", "MAT_Freq")], by = "Cell", all.x = T)
Ploidy_overview$Parent <- ifelse(Ploidy_overview$MAT_Freq < 0.15, "Paternal", "Mixed")
Ploidy_overview$Parent <- ifelse(Ploidy_overview$MAT_Freq > 0.5, "Maternal", Ploidy_overview$Parent)

## Ploidy based on Strand-seq
# Determine ploidy_SS based on distribution of strand-seq copy number states

StrandSeq_metrics$Haploidy <- (StrandSeq_metrics$somy0 + StrandSeq_metrics$somy1) / StrandSeq_metrics$number_segments

StrandSeq_metrics$Multiploidy <-  (StrandSeq_metrics$somy3 + StrandSeq_metrics$somy4) / StrandSeq_metrics$number_segments

StrandSeq_metrics[which(StrandSeq_metrics$WC_distribution < 0.2 & StrandSeq_metrics$Ploidy_SS == "Diploid"),]

StrandSeq_metrics$Ploidy_SS <- ifelse(StrandSeq_metrics$Haploidy > 0.9 & StrandSeq_metrics$WC_distribution < 0.25 & 
                                        StrandSeq_metrics$Method == "Strand-seq", 
                                         "Haploid", "ND")
StrandSeq_metrics$Ploidy_SS <- ifelse(StrandSeq_metrics$Haploidy > 0.5 & StrandSeq_metrics$Haploidy < 0.9 & 
                                        StrandSeq_metrics$WC_distribution < 0.4 & StrandSeq_metrics$Method == "Strand-seq", 
                                         "Largely_haploid", StrandSeq_metrics$Ploidy_SS)
StrandSeq_metrics$Ploidy_SS <- ifelse(StrandSeq_metrics$Haploidy > 0.3 & StrandSeq_metrics$Haploidy < 0.5 & 
                                        StrandSeq_metrics$WC_distribution < 0.4 & 
                                        StrandSeq_metrics$Method == "Strand-seq", 
                                         "Partially_haploid", StrandSeq_metrics$Ploidy_SS)
StrandSeq_metrics$Ploidy_SS <- ifelse(StrandSeq_metrics$Haploidy < 0.3 & 
                                        StrandSeq_metrics$Method == "Strand-seq", 
                                         "Diploid", StrandSeq_metrics$Ploidy_SS)

StrandSeq_metrics$Ploidy_SS <- ifelse(StrandSeq_metrics$Multiploidy > 0.75 & StrandSeq_metrics$Method == "Strand-seq", "Multiploid", StrandSeq_metrics$Ploidy_SS )

summary(factor(StrandSeq_metrics$Ploidy_SS))
# 
# SS_CN_state_overview$Library <- gsub(pattern = ".bam", x = SS_CN_state_overview$cell, replacement = "")
# SS_haploids <- gsub(pattern = ".bam", as.vector(SS_CN_state_overview$cell[which(SS_CN_state_overview$Ploidy == "Haploid")]), replacement = "")

Ploidy_overview <- merge(Ploidy_overview, StrandSeq_metrics[,c("Cell", "Ploidy_SS")], all.x = T)
Ploidy_overview$Ploidy_SS[is.na(Ploidy_overview$Ploidy_SS)] <- "ND"


## Ploidy conclusion
# Haploid/Uniparental if there are no/few heterozygous SNPs:
Ploidy_overview$Ploidy <- ifelse(Ploidy_overview$Ploidy_SNPs == "Haploid/Uniparental" & 
                                   Ploidy_overview$Ploidy_XY != "Diploid" & 
                                   Ploidy_overview$Ploidy_SS != "Diploid" & Ploidy_overview$Ploidy_SS != "Partially_haploid" & 
                                   Ploidy_overview$Parent != "Mixed" &
                                   Ploidy_overview$DEL0 %in% c("Partial", "ND"), "Haploid/Uniparental", "ND")

# Haploid if cell is Haploid/Uniparental based on het SNPs  AND deletions are nullisomies
Ploidy_overview$Ploidy <- ifelse(Ploidy_overview$Ploidy_SNPs == "Haploid/Uniparental" & 
                                   Ploidy_overview$Ploidy_XY != "Diploid" & 
                                   Ploidy_overview$Ploidy_SS != "Diploid" & Ploidy_overview$Ploidy_SS != "Partially_haploid" & 
                                   Ploidy_overview$DEL0 == "Yes", "Haploid", Ploidy_overview$Ploidy)

# Uniparental if there are no heterozygous SNPs, but deletions are not nullisomies (uniparental disomy)
Ploidy_overview$Ploidy <- ifelse(Ploidy_overview$Ploidy_SNPs == "Haploid/Uniparental" & 
                                   Ploidy_overview$Parent != "Mixed" &
                                   Ploidy_overview$Ploidy_XY != "Diploid" & 
                                   Ploidy_overview$Ploidy_SS != "Diploid" & Ploidy_overview$Ploidy_SS != "Partially_haploid" & 
                                   Ploidy_overview$DEL0 == "No" & Ploidy_overview$total.read.count > 100000, "Uniparental",Ploidy_overview$Ploidy)

# Uniparental if there are only maternally or only paternally inherited SNPs:
Ploidy_overview$Ploidy <- ifelse(Ploidy_overview$Parent != "Mixed" & 
                                   Ploidy_overview$Ploidy_SNPs != "Haploid/Uniparental" &
                                   Ploidy_overview$Ploidy_SS != "Diploid" & 
                                   Ploidy_overview$Ploidy_SS != "Partially_haploid" & 
                                   #Ploidy_overview$DEL0 != "Yes" & 
                                   Ploidy_overview$total.read.count > 50000, 
                                 "Uniparental",Ploidy_overview$Ploidy)

# 
Ploidy_overview$Ploidy <- ifelse(Ploidy_overview$Ploidy_SNPs != "Haploid/Uniparental" & 
                                   Ploidy_overview$Ploidy_SS == "Multiploid" & Ploidy_overview$DEL0 != "No" & Ploidy_overview$DEL0 != "Partial", "Multiploid", Ploidy_overview$Ploidy)

Ploidy_overview$Ploidy <- ifelse(Ploidy_overview$Ploidy_SS == "Haploid" & 
                                   Ploidy_overview$Ploidy_XY != "Diploid" , "Haploid", Ploidy_overview$Ploidy )

Ploidy_overview$Ploidy <- ifelse(Ploidy_overview$Ploidy_SS %in% c("Largely_haploid"), "Mixoploid", Ploidy_overview$Ploidy )


Ploidy_overview$Ploidy <- ifelse(Ploidy_overview$Ploidy_SNPs == "Diploid" & 
                                   Ploidy_overview$Ploidy_XY != "Haploid" & 
                                   Ploidy_overview$Ploidy_SS != "Haploid" & Ploidy_overview$Ploidy_SS != "Multiploid" &
                                   Ploidy_overview$Parent == "Mixed", "Diploid", Ploidy_overview$Ploidy )

summary(factor(Ploidy_overview$Ploidy ))


## Determine which cells are fragmented (fragmented cells have >25% of their reads on one chromosome)
Merged_idxstats_relative <- sweep(Merged_idxstats[,c(2:ncol(Merged_idxstats))],2,colSums(Merged_idxstats[,c(2:ncol(Merged_idxstats))]),`/`)
Merged_idxstats_relative[Merged_idxstats_relative < 0.25] <- 0
Fragmented_cells <- names(Merged_idxstats_relative[,which(colSums(Merged_idxstats_relative) > 0)])

Ploidy_overview$Ploidy[Ploidy_overview$Cell %in% Fragmented_cells] <- "Fragmented"

write.table(x = Ploidy_overview, file = paste(QC_folder, "Overview_Ploidy.txt", sep = ""), quote = F, sep = "\t")

# Plot SNPs etc
pdf(file = paste(results_folder, "Ploidy_analysis.pdf", sep = ""), width = 4, height = 3, pointsize = 10)

ggplot(Ploidy_overview[Ploidy_overview$total.read.count > 100000,], aes(x = AF0AF1, y = total.read.count, col = Ploidy_SNPs )) + geom_point(size = 0.5) + theme_bw(base_size = 8)

Ploidy_overview$Group <- factor(Ploidy_overview$Group)
ggplot(Ploidy_overview[Ploidy_overview$total.read.count > 150000,], aes(x = AF0AF1, fill = Ploidy_SNPs )) + geom_density(alpha = 0.7) + coord_cartesian(xlim = c(0,10)) + 
  labs(x = "Heterozygous SNPs / 1000 SNPs", fill = "Ploidy based on\nloss of heterozygosity") + theme_bw(base_size = 8)
ggplot(Ploidy_overview[Ploidy_overview$total.read.count > 150000,], aes(x = AF0AF1, fill = Ploidy_XY)) + geom_density(alpha = 0.7) + coord_cartesian(xlim = c(0,10)) +
  labs(x = "Heterozygous SNPs / 1000 SNPs", fill = "Ploidy based on\nabsence chrX") + theme_bw(base_size = 8)
ggplot(Ploidy_overview[which(Ploidy_overview$total.read.count > 150000 & Ploidy_overview$Method != "WGS"),], aes(x = AF0AF1, fill = Ploidy_SS)) + 
  geom_density(alpha = 0.6) + coord_cartesian(xlim = c(0,10)) + labs(x = "Heterozygous SNPs / 1000 SNPs", fill = "Ploidy based\non Strand-seq") + theme_bw(base_size = 8)

ggplot(Ploidy_overview[which(Ploidy_overview$total.read.count > 150000),], aes(fill = Ploidy, x = Group)) + geom_bar(stat = "count", position = "dodge")  + facet_grid(~ Stage) + theme_bw(base_size = 8)
  
dev.off()
