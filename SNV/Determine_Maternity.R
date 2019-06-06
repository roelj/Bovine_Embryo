# This script collects MAT and PAT/MAT SNP locations from a bed file and counts their frequencies per cell, chromosome or bin
library(GenomicRanges)
#library(ggplot2)
library(BSgenome.Btaurus.UCSC.bosTau8)
options(scipen = 999)

args <- commandArgs(trailingOnly=TRUE)

mode <- args[1] # cell | chr | bin
samplesheet_file <- args[2]
input_folder <- args[3] # requires bed files (chr - start - end - PAT/MAT)
output_folder <- args[4] # Subfolders will be created per run

input_parameters <- args[5] # This is used to set the bin size. separate by ";"
parameters <- unlist(strsplit(input_parameters, split = ";"))
samplesheet <- read.delim(samplesheet_file, stringsAsFactors = F)

# The SNP counts per bin will be written to this folder (one subfolder per run)
dir.create(output_folder, recursive = T)

cell_counts_overview <- data.frame()

# Length of the chromsomes is obtained from the BSgenome.Btaurus.UCSC.bosTau8 package
# This is required to count the maternal SP fraction per chromosome
# It is also used to generate genomic bins with a specified width
chr_sizes <- data.frame(chr = names(seqlengths(BSgenome.Btaurus.UCSC.bosTau8)[seqlevels(BSgenome.Btaurus.UCSC.bosTau8) %in% paste("chr", c(1:29, "X"), sep = "")]),
                        start = 1,
                        end = as.numeric(seqlengths(BSgenome.Btaurus.UCSC.bosTau8)[seqlevels(BSgenome.Btaurus.UCSC.bosTau8) %in% paste("chr", c(1:29, "X"), sep = "")]))
chr_sizes$chr <- gsub("chr", "", chr_sizes$chr)
chromosomes <- GRanges(seqnames = chr_sizes$chr, IRanges(start = chr_sizes$start, end = chr_sizes$end))

if(mode == "bin"){
  bin_size_input <-  parameters[grep(pattern = "bin_size", x = parameters)]
  if(length(bin_size_input) > 0 ){
    bin_size <- as.numeric(unlist(strsplit(bin_size_input, split = "="))[2])
  } else {
    print("# Bin_size not specified. Setting bin_size to 10000000.")
    bin_size <- 10e6
  } 
  # Generate genomic bins with specified width:
  bins <- tileGenome(seqlengths = seqlengths(BSgenome.Btaurus.UCSC.bosTau8)[seqlevels(BSgenome.Btaurus.UCSC.bosTau8) %in% paste("chr", c(1:29, "X"), sep = "")], tilewidth = bin_size, cut.last.tile.in.chrom=TRUE)
  seqlevels(bins) <-  gsub(seqlevels(bins), pattern = "chr", replacement = "")
}

# This function is used to calculate the maternal SNP fraction
# input_data should contain SNP locations with the following columns: chr - start - end - parent (MAT or PAT/MAT)
# Regions should be a GenomicRanges object or a string. If it is a string, the maternal fraction of the entire cell is calculated
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

print(paste("# Counting the number of maternal/paternal SNPs per ", mode, sep = ""))

for(Cell in samplesheet$Cell){
  #print(Cell)
  input_folder_run <- paste(input_folder, samplesheet$Run[samplesheet$Cell == Cell], "/", sep = "")
  input_file <- list.files(input_folder_run, full.names = T)[grep(pattern = paste(Cell,".bed", sep= ""), x = list.files(input_folder_run))]
  #print(input_file)
  if(length(input_file) > 0){
    if(file.size(input_file) > 0){
      input_data <- read.delim(input_file, stringsAsFactors = F, header = F)
      #print(head(input_data))
      if(mode == "cell"){
        freq_per_cell <- Count_SNPs_Bins(input_data = input_data, regions = "cell")
        cell_counts_overview <- rbind(cell_counts_overview, freq_per_cell)
        # output will be written below
      }
      
      if(mode == "chr"){
        # Calculate the maternal fraction per chromosome
        freq_per_chr <- Count_SNPs_Bins(input_data = input_data, regions = chromosomes)
        if(dir.exists(paste(output_folder, samplesheet$Run[samplesheet$Cell == Cell], "/", sep = "")) == FALSE){
          dir.create(paste(output_folder, samplesheet$Run[samplesheet$Cell == Cell], "/", sep = ""), recursive = T)
        }
        write.table(x = freq_per_chr, file =  paste(output_folder, samplesheet$Run[samplesheet$Cell == Cell], "/", Cell, ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
      }
      
      if(mode == "bin"){
        # Calculate the maternal fraction per 10 Mb bin
        freq_per_bin <- Count_SNPs_Bins(input_data = input_data, regions = bins)
        if(dir.exists(paste(output_folder, samplesheet$Run[samplesheet$Cell == Cell], "/", sep = "")) == FALSE){
          dir.create(paste(output_folder, samplesheet$Run[samplesheet$Cell == Cell], "/", sep = ""), recursive = T)
        }
        write.table(x = freq_per_bin, file = paste(output_folder, samplesheet$Run[samplesheet$Cell == Cell], "/", Cell, ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
      }
    }
  }
}

if(mode == "cell"){
  if(dir.exists(paste(output_folder, sep = "")) == FALSE){
    dir.create(paste(output_folder, sep = ""), recursive = T)
  }
  print(paste("# Writing: ", output_folder, "PAT_MAT_SNPs_per_cell.txt",sep = ""))
  write.table(cell_counts_overview, file = paste(output_folder, "PAT_MAT_SNPs_per_cell.txt", sep = ""), sep = "\t", quote = F, row.names = F)
}

print("# Done")
