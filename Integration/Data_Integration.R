# This scripts collects all the info about the cells and merges it into one overview
library(BSgenome.Btaurus.UCSC.bosTau8) #This package is only used to obtain the lengths of the chromosomes (which are used to determine aneuploidies)
library(GenomicRanges)
library(AneuFinder) # Only used here to load previously generated Aneufinder output

options(scipen = 999)
args <- commandArgs(trailingOnly=TRUE)
input_folder <- args[1] # [project_folder]

resolution <- 1000000

QC_folder <- paste(input_folder, "Data/QC/", sep = "")

samplesheet_file <- paste(QC_folder, "/Merged_quality_metrics.txt", sep = "")
samplesheet <-  read.delim(samplesheet_file, stringsAsFactors = F)

CNV_folder <-  paste(input_folder, "Data/CNV_Counts/", resolution, "/", sep = "")

# Add the number of CNVs to the sample overview
CNV_overview <- data.frame()
for(Embryo in unique(samplesheet$Embryo)){
  print(Embryo)
  CNV_file <- paste(CNV_folder, "CNV_Counts_", Embryo , ".txt", sep = "")
  
  if(file.exists(CNV_file)){
    CNVs_embryo <- read.delim(CNV_file, stringsAsFactors = F, check.names = F)
    CNV_counts <- data.frame(Cell = names(CNVs_embryo)[-c(1:2)], 
                             Aneuploidies = colSums(CNVs_embryo[which(CNVs_embryo$Type == "Aneuploidy" & CNVs_embryo$Ploidy %in% c("Gain", "Loss")),][,-c(1:2), drop=FALSE]),
                               CNVs = colSums(CNVs_embryo[which(CNVs_embryo$Type == "CNV" & CNVs_embryo$Ploidy %in% c("Gain", "Loss")),][,-c(1:2), drop=FALSE]),
                             Filtered_CNVs = colSums(CNVs_embryo[which(CNVs_embryo$Type == "Filter" & CNVs_embryo$Ploidy %in% c("Gain", "Loss", "Nullisomies")),][,-c(1:2), drop=FALSE]),
                             row.names = NULL)
    
    CNV_overview <- rbind(CNV_overview, CNV_counts)
  } else {
    CNV_counts <- data.frame(Cell = Cell, 
                             Aneuploidies = 0,
                             CNVs = 0,
                             Filtered_CNVs = 0)
    CNV_overview <- rbind(CNV_overview, CNV_counts)
    
  }
}

Cell_Overview <- merge(samplesheet, CNV_overview, all.x =T)
Cell_Overview$Filtered_CNVs[is.na(Cell_Overview$Filtered_CNVs)] <- 0
  
Ploidy_overview <- read.delim(paste(QC_folder, "Overview_Ploidy.txt", sep = ""), stringsAsFactors = F)
Cell_Overview <- merge(Cell_Overview , Ploidy_overview[,c("Cell", "Sex", "Parent", "Ploidy")], all.x = T)

## Cell conclusions
cell_overview <- Cell_Overview
cell_overview$Group <- factor(cell_overview$Group)
cell_overview$Conclusion <- ifelse(cell_overview$total.read.count > 100000, "Euploid", "ND")
cell_overview$Conclusion <- ifelse(cell_overview$CNVs > 0 & cell_overview$Conclusion == "Euploid", "CNV", cell_overview$Conclusion)
cell_overview$Conclusion <- ifelse(cell_overview$Aneuploidies > 0 & cell_overview$Conclusion %in% c("Euploid", "CNV"), "Aneuploid", cell_overview$Conclusion)
cell_overview$Conclusion <- ifelse(cell_overview$Ploidy == "Multiploid" & cell_overview$total.read.count > 100000, "Multiploid", cell_overview$Conclusion )
cell_overview$Conclusion <- ifelse(cell_overview$Ploidy == "Haploid" & !is.na(cell_overview$Ploidy), "Haploid", cell_overview$Conclusion )
cell_overview$Conclusion <- ifelse(cell_overview$Ploidy == "Haploid/Uniparental" & !is.na(cell_overview$Ploidy), "Haploid/Uniparental", cell_overview$Conclusion )
cell_overview$Conclusion <- ifelse(cell_overview$Ploidy == "Uniparental"& !is.na(cell_overview$Ploidy), "Uniparental", cell_overview$Conclusion )

cell_overview$Conclusion <- ifelse(cell_overview$Ploidy %in% c("Haploid", "Haploid/Uniparental", "Uniparental") & 
                                     !is.na(cell_overview$Ploidy) & cell_overview$CNVs + cell_overview$Aneuploidies > 2, 
                                   paste(cell_overview$Ploidy, "(Complex)"), 
                                   cell_overview$Conclusion)


cell_overview$Conclusion <- ifelse(cell_overview$Ploidy == "Fragmented" & !is.na(cell_overview$Ploidy) & cell_overview$total.read.count > 10000, "Fragmented", cell_overview$Conclusion )

cell_overview$Conclusion <- ifelse(cell_overview$CNVs + cell_overview$Aneuploidies > 2 & cell_overview$Conclusion %in% c("Euploid", "CNV", "Aneuploid"), "Complex", cell_overview$Conclusion)

summary(factor(cell_overview$base_cn_state))

cell_overview$base_cn_state[cell_overview$Conclusion %in% c("Haploid", "Haploid (Complex)") & cell_overview$base_cn_state < 3] <- 1

#cell_overview$Conclusion <- factor(cell_overview$Conclusion, levels = c("Complex", "Aneuploid","CNV", "Fragmented","Haploid","Haploid (Complex)", "Multiploid", "Euploid", "ND"))
summary(factor(cell_overview$Conclusion))

# cell exclusions (With reason). These cells are excluded from further analyses (for example counting of SVs)
cell_overview$Exclude <- "No"
cell_overview$Exclude[which(cell_overview$total.read.count < 100000)] <- "Low_read_count"
cell_overview$Exclude[which(cell_overview$Method == "WGS_Strand-seq" & 
                              cell_overview$Exclude != "No")] <- paste(cell_overview$Exclude[which(cell_overview$Method == "WGS_Strand-seq" & 
                                                                                                     cell_overview$Ploidy != "Fragmented" & 
                                                                                                     cell_overview$Exclude != "No")], "Failed_strand_seq", sep = ";")
cell_overview$Exclude[which(cell_overview$Method == "WGS_Strand-seq" & 
                              cell_overview$Exclude == "No")] <- "Failed_strand_seq"


cell_overview$Exclude[which(cell_overview$Filtered_CNVs > 10 & cell_overview$CNVs+cell_overview$Aneuploidies > 3 &
                              cell_overview$Exclude != "No")] <- paste(cell_overview$Exclude[which(cell_overview$Filtered_CNVs > 10 & cell_overview$CNVs+cell_overview$Aneuploidies > 3 &
                                                                                                     cell_overview$Exclude != "No")], "CNV_Filtered", sep = ";")
cell_overview$Exclude[which(cell_overview$Filtered_CNVs > 10 & cell_overview$CNVs+cell_overview$Aneuploidies > 3 &
                              cell_overview$Exclude == "No")] <- "CNV_Filtered"
cell_overview$Exclude[which(cell_overview$Conclusion == "Fragmented")] <- "No"



cell_overview$Exclude[which(cell_overview$num.segments > 80 &
                              cell_overview$Exclude != "No")] <- paste(cell_overview$Exclude[which(cell_overview$num.segments > 80 & 
                                                                                                     cell_overview$Exclude != "No")], "Many_segments", sep = ";")
cell_overview$Exclude[which(cell_overview$num.segments > 80 &
                              cell_overview$Exclude == "No")] <- "Many_segments"

cell_overview$Exclude[which(cell_overview$Conclusion == "Fragmented")] <- "No"



summary(factor(cell_overview$Exclude))

CNVs_embryo_folder <-  paste(input_folder, "Data/CNV/",resolution, "/", sep = "")


### Embryo classifications
print("## Classifying embryos")
print("# Determining number of different genotypes per embryo")
# Determine the number of genotypes per embryo
Embryo_overview <- data.frame()
for(embryo in unique(cell_overview$Embryo)){
  print(embryo)
  Embryo_conclusion <- data.frame(Embryo = embryo)
  Embryo_data <- cell_overview[cell_overview$Embryo == embryo,]
  conclusions <- Embryo_data$Conclusion[Embryo_data$Exclude == "No"]
  if(length(conclusions) > 1){
    #print(conclusions)
    Embryo_conclusion$Group <- unique(Embryo_data$Group)
    Embryo_conclusion$Stage <- unique(Embryo_data$Stage)
    Embryo_conclusion$Method <- ifelse(length(grep(pattern = "Strand-seq", x = Embryo_data$Method) > 0), "Strand-seq", "WGS")
    Embryo_conclusion$Run <- unique(Embryo_data$Run)
    
    Embryo_conclusion$Libraries <- nrow(Embryo_data)
    Embryo_conclusion$Cells <- length(which(Embryo_data$Exclude == "No"))
    Embryo_conclusion$Euploids <- length(which(conclusions == "Euploid"))
    
    if(file.exists(paste(CNVs_embryo_folder, "CNVs_", embryo, ".txt", sep = ""))){
      CNVs_Embryo <- read.delim(paste(CNVs_embryo_folder, "CNVs_", embryo, ".txt", sep = ""), check.names = F)
    
      Good_cells_embryo <- cell_overview$Cell[cell_overview$Embryo == embryo & cell_overview$Exclude == "No"]
    
    # Select all genotypes
    CNVs_Embryo <- CNVs_Embryo[which(CNVs_Embryo$width > 20e6),]
    if(nrow(CNVs_Embryo) == 0){
      Embryo_conclusion$Genotypes <-  1
    } else if(ncol(CNVs_Embryo) > 15){
      copy_numbers <- CNVs_Embryo[,which(names(CNVs_Embryo) %in% Good_cells_embryo)]
      copy_number_counts <- matrix(nrow = length(names(copy_numbers)), ncol = length(names(copy_numbers)), dimnames = list(names(copy_numbers), names(copy_numbers)))
      
      # Count the similarities and differences in CNVs between cells per pair of cells
      for(i in 1:(ncol(copy_numbers)-1)){
        #print(i)
        
        if(ncol(copy_numbers) > 2){
          for(n in i:ncol(copy_numbers[,2:ncol(copy_numbers)])){
            #print(n+1)
            pairs <- copy_numbers[,c(i,n+1)]
            if(length(which(pairs[,1] == 0 & pairs[,2] == 0)) > 0){
              pairs <- pairs[-which(pairs[,1] == 0 & pairs[,2] == 0),]
            }
            #print(names(pairs))
            
            similarities <- nrow(pairs[which(pairs[,1] < 0 & pairs[,2] < 0 | pairs[,1] > 0 & pairs[,2] > 0),])
            differences <- nrow(pairs) - similarities
            
            copy_number_counts[names(pairs)[1], names(pairs)[2]] <- similarities / (similarities+differences)
            copy_number_counts[names(pairs)[2], names(pairs)[1]] <- similarities / (similarities+differences)
          }
        } else {
          pairs <- copy_numbers
          #print(names(pairs))
          
          similarities <- nrow(pairs[which(pairs[,1] < 0 & pairs[,2] < 0 | pairs[,1] > 0 & pairs[,2] > 0),])
          differences <- nrow(pairs) - similarities
          
          copy_number_counts[names(pairs)[1], names(pairs)[2]] <- similarities / (similarities+differences)
          copy_number_counts[names(pairs)[2], names(pairs)[1]] <- similarities / (similarities+differences)
        }
    
      }
      
      copy_number_counts[is.nan(copy_number_counts)] <- 1
      groups2 <- ""
      
      # Group all cells based on their similarities. Cells are grouped if they share >75% of their CNVs
      for(cell in colnames(copy_number_counts)){
        #print(cell)
        cell_data <- copy_number_counts[cell,]
        group <- names(cell_data[cell_data >= 0.75])
        groups <- paste(c(cell, group[!is.na(group)])[order(c(cell, group[!is.na(group)]))], collapse = ",")
        groups2 <- c(groups2, groups)
      }
      groups2 <- groups2[groups2 != ""]
      Embryo_conclusion$Genotypes  <-  length(unique(groups2))
    }
    } else {
      Embryo_conclusion$Genotypes <- 1
    }
    Embryo_overview <- rbind(Embryo_overview, Embryo_conclusion) 
  }
}

Embryo_overview$Conclusion <- ifelse(Embryo_overview$Euploids/Embryo_overview$Cells >= 0.75, "Euploid", "Aneuploid/CNV")
Embryo_overview$Conclusion <- ifelse(Embryo_overview$Genotypes > 1 & Embryo_overview$Euploids/Embryo_overview$Cells < 0.75, "Mosaic", Embryo_overview$Conclusion)
Embryo_overview$Conclusion <- ifelse(Embryo_overview$Genotypes > 3, "Chaotic mosaic", Embryo_overview$Conclusion)
summary(factor(Embryo_overview$Conclusion))

Embryos_3cells <- as.character(Embryo_overview$Embryo[which(Embryo_overview$Stage == "2-cell" & Embryo_overview$Libraries > 2)])
Embryo_overview$Stage[Embryo_overview$Embryo %in% Embryos_3cells] <- "3-cell"
cell_overview$Stage[cell_overview$Embryo %in% Embryos_3cells] <- "3-cell"

for(embryo in cell_overview$Embryo){
  
  if(length( unique(Embryo_overview$Conclusion[Embryo_overview$Embryo == embryo])) > 0){
    cell_overview$Embryo_classification[cell_overview$Embryo == embryo] <- unique(Embryo_overview$Conclusion[Embryo_overview$Embryo == embryo])
    
  } else {
    cell_overview$Embryo_classification[cell_overview$Embryo == embryo] <- NA
  }
  
}

write.table(cell_overview, file = paste(QC_folder, "Overview_cells_",resolution,".txt", sep = ""), quote = F, sep = "\t", row.names = F)
write.table(Embryo_overview, file = paste(QC_folder, "Overview_embryos_",resolution,".txt", sep = ""), quote = F, sep = "\t", row.names = F)



## Determine copy number states per cell
print("## Determining the copy number states per genomic segment in each cell")
output_folder <- paste(input_folder, "Data/CN_States/", sep = "")
dir.create(output_folder)
export_figures <- FALSE

resolution <- 1000000

Overview_cells <- read.delim(paste(QC_folder, "Overview_cells_",resolution,".txt", sep = ""), stringsAsFactors = F)

# Read cnvs 
# The chromosome lengths are necessary to determine which % of a chromosome is affected by a CNV
chr_lengths <- seqlengths(BSgenome.Btaurus.UCSC.bosTau8)
names(chr_lengths) <- gsub(pattern = "chr", replacement = "", x = names(chr_lengths))
chromosomes <- c(1:29, "X", "Y")
chr_lengths <- chr_lengths[which(names(chr_lengths) %in% chromosomes)]
chr_lengths <- data.frame(seqnames= names(chr_lengths), chr_size = as.numeric(chr_lengths))
Genome <- GRanges(seqnames = chr_lengths$seqnames, IRanges(start = 1, end = chr_lengths$chr_size))

CNV_folder <- paste(input_folder, "Data/CNVs_cells/", resolution, "/", sep = "")
MAT_SNP_folder <-  paste(input_folder, "Data/SNVs/SNPs_Phased/", sep = "")


Count_SNPs_Bins <- function(input_data, regions){
  # MAT frequencies of the entire input data are calculated if "regions" is a string:
  if(is.character(regions) == TRUE){
    MAT_frequencies <- data.frame(Cell = Cell, MAT_Counts = length(which(input_data[,4] == "MAT")),
                                  PAT_MAT_Counts = length(which(input_data[,4] == "PAT/MAT")),
                                  MAT_Freq = length(which(input_data[,4] == "MAT"))/nrow(input_data),
                                  PAT_MAT_Freq = length(which(input_data[,4] == "PAT/MAT"))/nrow(input_data))
  } else {
    #print("regions")
    MAT_positions <- GRanges(seqnames = input_data[which(input_data[,4] == "MAT"),1], 
                             IRanges(start = input_data[which(input_data[,4] == "MAT"),2], end = input_data[which(input_data[,4] == "MAT"),3]))
    PAT_MAT_positions <- GRanges(seqnames = input_data[which(input_data[,4] == "PAT/MAT"),1], 
                                 IRanges(start = input_data[which(input_data[,4] == "PAT/MAT"),2], end = input_data[which(input_data[,4] == "PAT/MAT"),3]))
    
    MAT_Counts<- countOverlaps(regions, MAT_positions)
    PAT_MAT_Counts <- countOverlaps(regions, PAT_MAT_positions)
    
    MAT_frequencies <- data.frame(regions, MAT_Counts, PAT_MAT_Counts)
    MAT_frequencies$MAT_Freq <- MAT_frequencies$MAT_Counts / (MAT_frequencies$MAT_Counts+MAT_frequencies$PAT_MAT_Counts)
  }
  return(MAT_frequencies)
}

CNV_overview <- data.frame()
for(Cell in unique(Overview_cells$Cell)){
  print(Cell)
  Cell_overview <- Overview_cells[which(Overview_cells$Cell == Cell),]
  
  CN_states_phased <- data.frame(stringsAsFactors = F)
  CNV_file <- paste(CNV_folder, Cell_overview$Run, "/CNVs_", Cell , ".txt", sep = "")
  SNP_file <- paste(MAT_SNP_folder, Cell_overview$Run, "/", Cell , ".bed", sep = "")
  
  if(file.exists(CNV_file) == TRUE){
    CNVs_cell <- read.delim(CNV_file, stringsAsFactors = F, check.names = F)
    
    base <- ifelse(Cell_overview$base_cn_state > 1, Cell_overview$base_cn_state, 2)
    
    CNVs_cell$CN_state <- base + CNVs_cell[,Cell]
    if(Cell_overview$Ploidy == "Haploid"){
      CNVs_cell$CN_state <- ceiling(CNVs_cell$CN_state / 2)
    }
    
    CNVs_g <- GRanges(seqnames = CNVs_cell$seqnames, IRanges(start = CNVs_cell$start, end = CNVs_cell$end), 
                      Copy_number = CNVs_cell$CN_state)
  } else {
    CNVs_g <- GRanges()
  }
  
  Regions_base_CN <- setdiff(Genome, CNVs_g)
  Regions_base_CN$Copy_number <- Cell_overview$base_cn_state
  CN_cell <- c(CNVs_g,Regions_base_CN )
  if(file.exists(SNP_file) == TRUE){
    if(file.size(SNP_file) > 0){
      SNPs_cell <- read.delim(SNP_file, stringsAsFactors = F, check.names = F)
      CN_states_phased <- Count_SNPs_Bins(input_data = SNPs_cell, regions = CN_cell)
    } else {
      CN_states_phased <- data.frame(CN_cell)
      CN_states_phased$MAT_Counts <- NA
      CN_states_phased$PAT_MAT_Counts  <- NA
      CN_states_phased$MAT_Freq <- NA
    }
  }  else {
    CN_states_phased <- data.frame(CN_cell)
    CN_states_phased$MAT_Counts <- NA
    CN_states_phased$PAT_MAT_Counts  <- NA
    CN_states_phased$MAT_Freq <- NA
  }
  CN_states_phased <- merge(CN_states_phased, chr_lengths, by = "seqnames", all.x =T)
  CN_states_phased$Relative_width <- round(CN_states_phased$width / CN_states_phased$chr_size, 2)
  output_folder_run <- paste(output_folder, unique(Cell_overview$Run), "/", sep = "")
  if(dir.exists(output_folder_run) == FALSE){
    dir.create(output_folder_run, recursive = T)
  }
  write.table(CN_states_phased, file = paste(output_folder_run, "CN_States_", Cell, ".txt", sep = ""), quote = F, sep = "\t", row.names = F)
}


