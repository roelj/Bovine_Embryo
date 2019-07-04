## This script integrates copy number variant calls of each cell with an embryo.
# CNV calls are merged based on a certain percentage of overlap depending on the size of the CNV.
# One callset per embryo is generated and subsequently the presence of the merged CNVs in the single cells is determined.

#source("https://bioconductor.org/biocLite.R")
#BiocInstaller::biocLite(c("BSgenome.Btaurus.UCSC.bosTau8"))
library(BSgenome.Btaurus.UCSC.bosTau8) #This package is only used to obtain the lengths of the chromosomes (which are used to determine aneuploidies)
library(GenomicRanges)
library(AneuFinder) # Only used here to load previously generated Aneufinder output
options(scipen = 999)

args <- commandArgs(trailingOnly=TRUE)
input_folder <- args[1]
resolution <- args[2]
blacklist_path <- args[3]
#output_folder <- args[2]
runs <- args[4]

# Settings
QC_folder <- paste(input_folder, "Data/QC/", sep = "")

samplesheet_file <- paste(QC_folder, "/Merged_quality_metrics.txt", sep = "") # Output of Collect_QualityMetrics.R script
aneufinder_folder <-  paste(input_folder,"Data/Aneufinder/", sep = "") # This folder should contain the output of Aneufinder (one output folder per run)
output_folder <-  paste(input_folder, "Data/CNV/", resolution,"/", sep = "") # Filtered CNVs will be written here
CNV_Counts_Folder <- paste(input_folder, "Data/CNV_Counts/", resolution, "/", sep = "") # CNV counts per embryo will be written here
dir.create(CNV_Counts_Folder, showWarnings = F, recursive = T)
#blacklist_path <- "" # Leave empty if you want to generate a new blacklist
#blacklist_path <- paste(output_folder, "CNV_blacklist.txt", sep = "")
#runs <- c("") # Leave empty if you want to run for all runs/folders
make_plots <- FALSE # Set to true to plot population frequencies of CNVs

embryo_overview <- data.frame()
# Load settings
print(paste("# Reading samplesheet: ", samplesheet_file, sep = ""))
samplesheet <- read.delim(samplesheet_file, header = T, stringsAsFactors = F)
# Remove positive and negative controls:
samplesheet <- samplesheet[samplesheet$Embryo != "POS" & samplesheet$Embryo != "NEG",]

resolution <- as.numeric(resolution)

# Only select the runs defined in runs:
if(runs != ""){
  runs <- unlist(strsplit(runs, split = ","))
  samplesheet <- samplesheet[samplesheet$Run %in% runs,]
}
dir.create(output_folder, recursive = T)

# The chromosome lengths are necessary to determine which % of a chromosome is affected by a CNV
chr_lengths <- seqlengths(BSgenome.Btaurus.UCSC.bosTau8)
names(chr_lengths) <- gsub(pattern = "chr", replacement = "", x = names(chr_lengths))
chromosomes <- c(1:29, "X", "Y")
chr_lengths <- chr_lengths[which(names(chr_lengths) %in% chromosomes)]
chr_lengths <- data.frame(seqnames= names(chr_lengths), chr_size = as.numeric(chr_lengths))

## Define functions
# This function can be used to merge regions (CNVs) based on a percentage of overlap between the regions (x and y, Genomicranges object). 
# This percentage is adjusted for the size of the region. For example, it increases for larger regions.
# The function can also be used to filter away overlapping regions, by setting ' mode = "Filter" '
Merge_Overlapping_CNVs <- function(x, y, minimal_overlap = 0.3, maximum_overlap = 0.9, mode = "Merge") {
  x <- granges(x)
  y <- granges(y)
  hits <- findOverlaps(x, y)
  
  xhits <- x[queryHits(hits)]
  yhits <- y[subjectHits(hits)]
  
  # determine how much overlap there is between the two ranges:
  frac <- width(pintersect(xhits, yhits)) / pmax(width(xhits), width(yhits))
  # determine the minimal overlap that is necessary to merge the two ranges.
  # depends on the size of the largest fragment
  minfrac <- minimal_overlap + floor(pmax(width(xhits), width(yhits)) / 10e6) * 0.1
  minfrac[which(minfrac > maximum_overlap)] <- maximum_overlap
  
  merge <- frac >= minfrac
  
  if(mode == "Merge"){
      # The CNVs that overlap with more than the minfrac will be merged together.
      merged <- c(reduce(c(xhits[merge], yhits[merge])),
                  xhits[!merge], yhits[!merge],
                  x[-queryHits(hits)], y[-subjectHits(hits)])
      # The merging generates duplicates, remove the duplicates:
      merged <- merged[!duplicated(merged)]
      } else if (mode == "Filter"){
        if(length(hits) > 0){
      # Remove the entries in x that overlap with y
      merged <-  c(x[-queryHits(hits)], xhits[!merge])
      } else {
        merged <- x
    }
  }

  return(merged)
}

# This function counts the number of overlapping regions between x and y. 
# The required overlap is adjusted for region size (the larger the region the higher percentage of overlap is required)
Count_Overlapping_CNVs <- function(x, y, minimal_overlap=0.3, maximum_overlap = 0.9) {
  #x <- granges(x)
  #y <- granges(y)
  hits <- findOverlaps(x, y)
  xhits <- x[queryHits(hits)]
  yhits <- y[subjectHits(hits)]
  frac <- width(pintersect(xhits, yhits)) / pmax(width(xhits), width(yhits))
  
  minfrac <- minimal_overlap + floor(pmax(width(xhits), width(yhits)) / 10e6) * 0.1
  minfrac[which(minfrac > maximum_overlap)] <- maximum_overlap
  
  merge <- frac >= minfrac
  
  xhits$counts <- yhits$counts
  #xhits$counts[merge] <- 1
  
  return(xhits[merge])
}

# This function first selects CNVs with the same copy number state and subsequently extends these CNVs with a certain window.
# Subsequently the extended CNVs are reduces (collapsed) into larger CNVs (if there is overlap between extended CNVs)
# Next the raw, unextended CNVs are again overlapped with the extended CNVs.
# The minimal and maximal start and end positions of CNVs overlapping with extended CNVs are used to determine the start and end of chained CNVs.
Chain_CNVs <- function(Input_GRanges, Window_percentage = 0.2, Max_window = 7.5e6){
  output <- GRanges()
  
  # loop over each copy number state in the input set to only merge CNVs of same copy number:
  for(cn_state in unique(Input_GRanges$copy.number)){
    #print(cn_state)
    
    # select the CNVs with the same cn state:
    CNVs_state <- Input_GRanges[which(Input_GRanges$copy.number == cn_state),]
    
    Extended_CNVs <- CNVs_state
    #print(CNVs_state)
    
    # Extend the selected CNVs with 10% (Window_percentage) of their width (max 10 Mb, Max_window)
    for(i in 1:length(CNVs_state)){
      window <- width(Extended_CNVs)[i]*Window_percentage
      window <- ifelse(window > Max_window, Max_window, window)
      start(Extended_CNVs)[i] <- ifelse((start(Extended_CNVs)[i] - window < 1), 1, (start(Extended_CNVs)[i] - window))
      end(Extended_CNVs)[i] <- end(Extended_CNVs)[i] + window
    }
    Extended_CNVs <- trim(Extended_CNVs)
    
    # Reduce the extended CNVs to create large CNVs
    Reduced_extended_CNVs <- reduce(Extended_CNVs)
    
    # Overlap the raw CNVs with the extended large CNVs:
    Overlap_raw_extended_CNVs <- findOverlaps(CNVs_state, Reduced_extended_CNVs)
    
    # Remove the extensions up and downstream of the chained CNVs:
    for(overlap in unique(subjectHits(Overlap_raw_extended_CNVs))){
      overlapping_regions <- CNVs_state[queryHits(Overlap_raw_extended_CNVs)[subjectHits(Overlap_raw_extended_CNVs) == overlap]]
      
      output_overlap <- GRanges(seqnames = unique(seqnames(overlapping_regions)), 
                                IRanges(start = min(start(overlapping_regions)), end = max(end(overlapping_regions)) ),
                                copy.number = cn_state)
      output <- c(output, output_overlap)
      output <- output[!duplicated(output)]
    }
  }
  
  # Remove overlapping events
  output$ID <- paste(seqnames(output), start(output), end(output), sep = "_")
  output2 <- output
  # Add 1 to the start positions to prevent merging of adjacent regions by reduce()
  start(output2) <- start(output2)+1
  output2 <- reduce(output2)
  start(output2) <- start(output2)-1
  output2$ID <- paste(seqnames(output2), start(output2), end(output2), sep = "_")
  
  # Only keep the non-overlapping events (largest is kept)
  output <- output[output$ID %in% output2$ID]
  return(output)
}

# Plot the threshold:
if(make_plots == TRUE){
  threshold_overview <- data.frame(Size = seq(0, 120, by = 1))
  threshold_overview$Percentage <- threshold_overview$Size / 10 * 0.1 * 100
  threshold_overview$Percentage <- ifelse(threshold_overview$Percentage < 30, 30, threshold_overview$Percentage)
  threshold_overview$Percentage <- ifelse(threshold_overview$Percentage > 90, 90, threshold_overview$Percentage)
  
  ggplot(threshold_overview, aes(x = Size, y = Percentage)) + geom_line(lwd = 2) + 
    scale_y_continuous(limits = c(0,100), expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) + labs(x = "Size region (Mb)", y = "Overlap threshold (%)") +
    theme(panel.grid.major = element_line(colour = "lightgray", linetype = 2))
  
  ggplot(threshold_overview, aes(x = Size, y = Percentage/100*Size)) + geom_line(lwd = 2) + 
    scale_y_continuous(limits = c(0,125), expand = c(0,0)) + 
    scale_x_continuous(limits = c(0,125), expand = c(0,0)) + labs(x = "Size region (Mb)", y = "Required overlap (Mb)") +
    theme(panel.grid.major = element_line(colour = "lightgray", linetype = 2))
}

## Blacklist
# First merge the CNVs of "normal" untreated embryos to generate a blacklist of "population" CNVs
if(blacklist_path == ""){
  samples_blacklist <- samplesheet[which(samplesheet$Method == "WGS" & samplesheet$Group == 0 & samplesheet$total.read.count > 200000),]
  merged_CNVs_population <- GRanges()
  merged_CNVs_population_post <- GRanges()
  
  print(paste("## Generating CNV blacklist based on ", nrow(samples_blacklist), " samples from ", length(unique(samples_blacklist$Run)), " runs", sep = ""))
  for(run in unique(samples_blacklist$Run)){
    print(run)
    aneufinder_output <- ifelse(length(grep("Strand-seq", x = unique(samplesheet[samplesheet$Run == run, "Method"])) > 0),
                                paste(aneufinder_folder, run, "/MODELS-StrandSeq/method-edivisive/", sep = ""),
                                paste(aneufinder_folder, run, "/MODELS/method-edivisive/", sep = ""))
    aneufinder_files <- list.files(aneufinder_output, full.names = T)
    print(paste("## Loading ", length(aneufinder_files[grep(pattern = paste("bam_binsize_", format(resolution, scientific = T), sep = ""), x = as.character(aneufinder_files), fixed = T)]), " Aneufinder files from folder: ", aneufinder_output, sep = ""))
    
    aneufinder_data <- loadFromFiles(aneufinder_files[grep(pattern = paste("bam_binsize_", format(resolution, scientific = T), sep = ""), x = as.character(aneufinder_files), fixed = T)])
    
    for(sample in samples_blacklist$Cell[samples_blacklist$Run == run]){
      print(sample)
      sample_info <- samplesheet[which(samplesheet$Cell == sample),]
      
      if(sample_info$total.read.count > 150000){
        
        aneufinder_data_cell <- aneufinder_data[[names(aneufinder_data)[grep(paste(sample, ".bam", sep = ""), x = names(aneufinder_data))]]]$segments
        
        if(!is.null(aneufinder_data_cell)){
          
          cnvs_cell <- aneufinder_data_cell[which(aneufinder_data_cell$copy.number != sample_info$base_cn_state),]
          
          if(length(merged_CNVs_population) == 0){
            merged_CNVs_population <- cnvs_cell
          } else if (length(merged_CNVs_population_post) == 0) {
            merged_CNVs_population <- Merge_Overlapping_CNVs(merged_CNVs_population, cnvs_cell, minimal_overlap = 0.3)
          } else {
            merged_CNVs_population <- Merge_Overlapping_CNVs(merged_CNVs_population_post, cnvs_cell, minimal_overlap = 0.3)
            
          }
        }
      }
    }
    merged_CNVs_population_post <- merged_CNVs_population
    print("Post filtering merged CNVs")
    
    for(i in length(merged_CNVs_population):1){
      merged_CNVs_population_post <- Merge_Overlapping_CNVs(merged_CNVs_population_post, merged_CNVs_population[i], minimal_overlap = 0.4)
    }
    
    postfiltered <- length(merged_CNVs_population) - length(merged_CNVs_population_post)
    
    print(paste("Post filtering removed ", postfiltered, " of " ,length(merged_CNVs_population), " merged CNVs" ,sep = ""))
  }
  
  i <- 0
  overview_population <- data.frame(merged_CNVs_population_post)
  overview_population$CNV_ID <- paste(overview_population$seqnames, overview_population$start, overview_population$end, sep = "_")
  for(run in unique(samples_blacklist$Run)){
    print(run)
    aneufinder_output <- ifelse(length(grep("Strand-seq", x = unique(samplesheet[samplesheet$Run == run, "Method"])) > 0),
                                paste(aneufinder_folder, run, "/MODELS-StrandSeq/method-edivisive/", sep = ""),
                                paste(aneufinder_folder, run, "/MODELS/method-edivisive/", sep = ""))
    aneufinder_files <- list.files(aneufinder_output, full.names = T)
    print(paste("## Loading ", length(aneufinder_files[grep(pattern = paste("bam_binsize_", format(resolution, scientific = T), sep = ""), x = as.character(aneufinder_files), fixed = T)]), " Aneufinder files from folder: ", aneufinder_output, sep = ""))
    
    aneufinder_data <- loadFromFiles(aneufinder_files[grep(pattern = paste("bam_binsize_", format(resolution, scientific = T), sep = ""), x = as.character(aneufinder_files), fixed = T)])
    for(sample in samples_blacklist$Cell[samples_blacklist$Run == run]){
      print(sample)
      sample_info <- samplesheet[which(samplesheet$Cell == sample),]
      
      if(sample_info$total.read.count > 150000){
        
        aneufinder_data_cell <- aneufinder_data[[names(aneufinder_data)[grep(paste(sample, ".bam", sep = ""), x = names(aneufinder_data))]]]$segments
        
        if(!is.null(aneufinder_data_cell)){
          cnvs_cell <- aneufinder_data_cell[which(aneufinder_data_cell$copy.number != 2),]
          
          cnvs_cell$counts <- ifelse(cnvs_cell$copy.number > sample_info$base_cn_state, 1, -1)
          cnvs_cell$counts <- ifelse(cnvs_cell$copy.number > sample_info$base_cn_state + 1, 2, cnvs_cell$counts)
          cnvs_cell$counts <- ifelse(cnvs_cell$copy.number < sample_info$base_cn_state -1, -2, cnvs_cell$counts)
          
          cnv_count_cell <- data.frame(Count_Overlapping_CNVs(x = merged_CNVs_population_post, y = cnvs_cell, minimal_overlap = 0.3))
          cnv_count_cell$CNV_ID <- paste(cnv_count_cell$seqnames, cnv_count_cell$start, cnv_count_cell$end, sep = "_")
          
          overview_population <- merge(overview_population, cnv_count_cell[,c("CNV_ID", "counts")], all.x = T, by = "CNV_ID")
          names(overview_population)[ncol(overview_population)] <- sample
          overview_population <- overview_population[order(overview_population$seqnames, overview_population$start),]
          i <- i + 1
          
        }
      }
    }
  }
  overview_population[is.na(overview_population)] <- 0
  overview_population$Population_count <- 0
  for(i in 1:nrow(overview_population)){
    
    overview_population$Population_count[i] <- length(which(overview_population[i,-c(1:6)] != 0))
  }
  
  overview_population$Pop_freq <- round( overview_population$Population_count / (ncol(overview_population) -  6),2)
  blacklist <- overview_population[which(overview_population$Pop_freq > 0.15 & overview_population$seqnames != "X"),]
  
  write.table(blacklist, file = paste(output_folder, "CNV_blacklist.txt", sep = ""), sep = "\t", quote = F, row.names = F)
  
} else {
  blacklist <- read.delim(blacklist_path, stringsAsFactors = F)
}

if(make_plots == TRUE){
  ggplot(overview_population, aes(x = Pop_freq)) + geom_density(fill = "gray") + scale_y_continuous(expand = c(0,0)) + labs(x = "Population frequency CNVs")
}



## Determine dynamic threshold
# Here we determine in which cell division an CNV or aneuploidy occurred. 
# Because not all embryos have the same number of cell sequenced we need a dynamic threshold for each cell division.
threshold <- list()

for(number_cells in 4:8){
  thresholds_division <- data.frame(Division = c(0:3, "2-3", "MN"))
  
  # Reciprocal and non-reciprocal CNVs have different thresholds.
  thresholds_division$Non_reciprocal_min <- round(ifelse(number_cells > 6, ceiling(number_cells-1)/number_cells, 1), 2)
  thresholds_division$Non_reciprocal_min[2] <- round(ifelse(number_cells > 4, ceiling(number_cells/2-1)/number_cells, number_cells/2/number_cells), 2)
  thresholds_division$Non_reciprocal_min[3] <- round(ifelse(number_cells > 5, ceiling(number_cells/4)/number_cells, NA), 2)
  thresholds_division$Non_reciprocal_min[4] <- round(ifelse(number_cells > 5, 1/number_cells, NA), 2)
  thresholds_division$Non_reciprocal_min[5] <- round(ifelse(number_cells > 5, NA, 1/number_cells), 2)
  
  thresholds_division$Non_reciprocal_max <- 1
  thresholds_division$Non_reciprocal_max[2] <-  round(ifelse(number_cells > 4,ceiling(number_cells/2)/number_cells, 0), 2)
  thresholds_division$Non_reciprocal_max[3] <-  round(ifelse(number_cells > 5,ceiling(number_cells/4)/number_cells, NA), 2)
  thresholds_division$Non_reciprocal_max[4] <-  round(ifelse(number_cells > 5,ceiling(number_cells/4)/number_cells, NA), 2)
  thresholds_division$Non_reciprocal_max[5] <- round(ifelse(number_cells > 5, NA, 1/number_cells), 2)
  
  # seperate category for events that occur 5/6 or 5 times in 8, 7, 6 cell embryos
  thresholds_division$Non_reciprocal_min[6] <- round(ifelse(number_cells > 5, thresholds_division$Non_reciprocal_max[2]+0.01, NA), 2)
  
  thresholds_division$Non_reciprocal_max[6] <- round(ifelse(number_cells > 5, thresholds_division$Non_reciprocal_min[1], NA), 2)
  
  
  ## The thresholds are different for reciprocal events.
  # for example: by definition a reciprocal event should always occur twice or more
  thresholds_division$Reciprocal_min <- NA
  thresholds_division$Reciprocal_min[2] <- round(ifelse(number_cells > 4, ceiling(number_cells-1)/number_cells, number_cells), 2)
  
  thresholds_division$Reciprocal_min[3] <- round(ifelse(number_cells > 4, ceiling(number_cells/2-1)/number_cells, number_cells/2/number_cells), 2)
  thresholds_division$Reciprocal_min[4] <- round(ifelse(number_cells > 5, ceiling(number_cells/4)/number_cells, NA), 2)
  thresholds_division$Reciprocal_min[5] <- NA
  
  thresholds_division$Reciprocal_max <- NA
  thresholds_division$Reciprocal_max[2] <- 1
  thresholds_division$Reciprocal_max[3] <-  round(ifelse(number_cells > 4,ceiling(number_cells/2)/number_cells, 0), 2)
  thresholds_division$Reciprocal_max[4] <-  round(ifelse(number_cells > 5,ceiling(number_cells/4)/number_cells, NA), 2)
  thresholds_division$Reciprocal_max[5] <- NA
  
  thresholds_division$Reciprocal_min[6] <- round(ifelse(number_cells > 6, thresholds_division$Reciprocal_max[3]+0.01, NA), 2)
  thresholds_division$Reciprocal_max[6] <- round(ifelse(number_cells > 6, thresholds_division$Reciprocal_min[2], NA), 2)
  
  threshold[[number_cells]] <- thresholds_division
}



## Step 1
# Merge the CNV calls for each embryo into one call set per embryo. Count the CNV calls per cell in the embryo
for(run in unique(samplesheet$Run)){
    
  print(run)
  aneufinder_output <- ifelse(length(grep("Strand-seq", x = unique(samplesheet[samplesheet$Run == run, "Method"])) > 0),
                              paste(aneufinder_folder, run, "/MODELS-StrandSeq/method-edivisive/", sep = ""),
                              paste(aneufinder_folder, run, "/MODELS/method-edivisive/", sep = ""))
  
  aneufinder_files <- list.files(aneufinder_output, full.names = T)
  
  print(paste("## Loading ", length(aneufinder_files[grep(pattern = paste("bam_binsize_", format(resolution, scientific = T), sep = ""), x = as.character(aneufinder_files), fixed = T)]), " Aneufinder files from folder: ", aneufinder_output, sep = ""))
  
  aneufinder_data <- loadFromFiles(aneufinder_files[grep(pattern = paste("bam_binsize_", format(resolution, scientific = T), sep = ""), x = as.character(aneufinder_files), fixed = T)])
  
    
    for(embryo in unique(samplesheet$Embryo[samplesheet$Run == run])){
      
      # The number of CNVs that do not pass the filter are counted for each cell and these counts are collected in the filtered_events list
      filtered_events <- list()
      
      print(paste("## Merging CNVs for embryo: ",embryo ,sep = ""))
      merged_CNVs_embryo <-  GRanges()
      embryo_overview <- data.frame()
      CNVs_filter_pass <- list()
      for(sample in samplesheet$Cell[samplesheet$Embryo == embryo]){
        #print(sample)
        sample_info <- samplesheet[which(samplesheet$Cell == sample),]
        
        # Only include samples with > 100.000 reads
        if(sample_info$total.read.count > 100000){
          
          # Variants in aneufinder_data_cell will be used for merging and counting. 
          # Variants in aneufinder_data_cell_filtered will be more strictly filtered and used to filter out unique SVs after merging and counting
          aneufinder_data_cell <- aneufinder_data[[names(aneufinder_data)[grep(paste(sample, ".bam", sep = ""), x = names(aneufinder_data))]]]$segments
          CNVs_filter_pass_cell <- aneufinder_data_cell
          
          # filter cnvs that have too little difference with the copy number state 2
          counts_2state <- mean(aneufinder_data_cell$mean.counts[aneufinder_data_cell$copy.number == sample_info$base_cn_state])
          
          if(length(aneufinder_data_cell$copy.number[which(aneufinder_data_cell$copy.number == 0 & aneufinder_data_cell$mean.counts > counts_2state/4)]) >0){
            filtered_events[[sample]] <- data.frame(ID = sample, 
                                                    Filtered_nullisomies = length(aneufinder_data_cell$copy.number[which(aneufinder_data_cell$copy.number == 0 & aneufinder_data_cell$mean.counts > counts_2state/4)]))
          } else {
            filtered_events[[sample]] <- data.frame(ID = sample, Filtered_nullisomies = 0)
          }

          # Set false nullisomies to copy state 1
          aneufinder_data_cell$copy.number[which(aneufinder_data_cell$copy.number == 0 & aneufinder_data_cell$mean.counts > counts_2state/4)] <- 1
          aneufinder_data_cell$state[which(aneufinder_data_cell$copy.number == 0 & aneufinder_data_cell$mean.counts > counts_2state/4)] <- "1-somy"
          
           # Recurrent CNVs within one embryo must have less or more reads than these filters:
          loss_threshold <- 0.75
          gain_threshold <- 1.25
          
          # Unique CNVs in an embryo will be more stringently filtered to reduce amount of false positives. 
          # These will be filtered after merging and counting the CNVs
          loss_threshold_stringent <- ifelse(sample_info$base_cn_state < 3, 0.65, 0.75) # loss from 3 to 2 is just a 33% change
          gain_threshold_stringent <- ifelse(sample_info$base_cn_state < 3, 1.4, 1.2) # gain from 3 to 4 is just a 25% change
          
          # Remove the deletions that do not pass the loss_threshold (eg mean counts in more than counts_2state*loss_threshold)
          if(length(which(aneufinder_data_cell$copy.number < sample_info$base_cn_state & aneufinder_data_cell$mean.counts > counts_2state*loss_threshold)) > 0){
            # filtered_events[[sample]] <- data.frame(ID = sample, Filtered_losses = length(aneufinder_data_cell[which(aneufinder_data_cell$copy.number < sample_info$base_cn_state & 
            #                                                                                                             aneufinder_data_cell$mean.counts > counts_2state*loss_threshold)]))
            aneufinder_data_cell <- aneufinder_data_cell[-which(aneufinder_data_cell$copy.number < sample_info$base_cn_state & aneufinder_data_cell$mean.counts > counts_2state*loss_threshold)]
          } 
          if(length(which(aneufinder_data_cell$copy.number > sample_info$base_cn_state & aneufinder_data_cell$mean.counts < counts_2state*gain_threshold)) > 0){
            #filtered_events[[sample]]$Filtered_gains <- length(aneufinder_data_cell[which(aneufinder_data_cell$copy.number > sample_info$base_cn_state & aneufinder_data_cell$mean.counts < counts_2state*gain_threshold)])
            aneufinder_data_cell <- aneufinder_data_cell[-which(aneufinder_data_cell$copy.number > sample_info$base_cn_state & aneufinder_data_cell$mean.counts < counts_2state*gain_threshold)]
          } 
          
          # The CNVs that do not pass the stringent filter will be stored in the CNVs_filter_pass list
          #CNVs_filter_pass_cell <- Chain_CNVs(CNVs_filter_pass_cell)
          
          CNVs_filter_pass_cell <- CNVs_filter_pass_cell[which(CNVs_filter_pass_cell$copy.number != sample_info$base_cn_state)]
          
          if(length(CNVs_filter_pass_cell) > 0){
            if(length(which(CNVs_filter_pass_cell$copy.number < sample_info$base_cn_state & CNVs_filter_pass_cell$mean.counts > counts_2state*loss_threshold_stringent) > 0)){
              
              filtered_events[[sample]]$Filtered_losses <- length(CNVs_filter_pass_cell[which(CNVs_filter_pass_cell$copy.number < sample_info$base_cn_state & 
                                                                                                CNVs_filter_pass_cell$mean.counts > counts_2state*loss_threshold_stringent)])
              
              CNVs_filter_pass_cell <- CNVs_filter_pass_cell[-which(CNVs_filter_pass_cell$copy.number < sample_info$base_cn_state & CNVs_filter_pass_cell$mean.counts > counts_2state*loss_threshold_stringent)]
                 
            } else {
              filtered_events[[sample]]$Filtered_losses <- 0
            }
            
            if(length(which(CNVs_filter_pass_cell$copy.number > sample_info$base_cn_state & CNVs_filter_pass_cell$mean.counts < counts_2state*gain_threshold_stringent) > 0)){
              filtered_events[[sample]]$Filtered_gains <- length(CNVs_filter_pass_cell[which(CNVs_filter_pass_cell$copy.number > sample_info$base_cn_state &CNVs_filter_pass_cell$mean.counts < counts_2state*gain_threshold_stringent)])
              
               CNVs_filter_pass_cell <- CNVs_filter_pass_cell[-which(CNVs_filter_pass_cell$copy.number > sample_info$base_cn_state & CNVs_filter_pass_cell$mean.counts < counts_2state*gain_threshold_stringent)]
            } else {
              filtered_events[[sample]]$Filtered_gains <- 0
            }
          }
          
          CNVs_filter_pass_cell <- Chain_CNVs(CNVs_filter_pass_cell)
          
          CNVs_filter_pass[[sample]] <- CNVs_filter_pass_cell
          
          if(!is.null(aneufinder_data_cell)){
            
            # Chain adjacent CNVs
            cnvs_cell <- aneufinder_data_cell[which(aneufinder_data_cell$copy.number != sample_info$base_cn_state),]
            cnvs_cell <- Chain_CNVs(cnvs_cell)
            end(cnvs_cell) <- end(cnvs_cell) - 1 
            #print(data.frame(cnvs_cell))
            if(length(merged_CNVs_embryo) == 0){
              merged_CNVs_embryo <- cnvs_cell
            } else {
              merged_CNVs_embryo <- Merge_Overlapping_CNVs(merged_CNVs_embryo, cnvs_cell, minimal_overlap = 0.5, maximum_overlap = 0.9)
              #print(data.frame(merged_CNVs_embryo))
            }
          }
        }
      }
      end(merged_CNVs_embryo) <- end(merged_CNVs_embryo)+1
      
      
      merged_CNVs_embryo_post <- merged_CNVs_embryo
      if(length(merged_CNVs_embryo_post) > 0){
        print(paste("Post filtering merged CNVs embryo: ",embryo ,sep = ""))
        
        for(i in length(merged_CNVs_embryo):1){
          merged_CNVs_embryo_post <- Merge_Overlapping_CNVs(merged_CNVs_embryo_post, merged_CNVs_embryo[i], minimal_overlap = 0.5, maximum_overlap = 0.9)
        }
        
        ## Remove CNVs that frequently occur in the population (reference errors / polymorphisms)
        merged_CNVs_embryo_post <- Merge_Overlapping_CNVs(x = merged_CNVs_embryo_post, 
                                                          y = GRanges(blacklist$seqnames, IRanges(start = blacklist$start, end = blacklist$end)), 
                                                          mode = "Filter", minimal_overlap = 0.5, maximum_overlap = 0.9)
        
        
        postfiltered <- length(merged_CNVs_embryo) - length(merged_CNVs_embryo_post)
        
        print(paste("Post filtering removed ", postfiltered, " of " ,length(merged_CNVs_embryo), " merged CNVs" ,sep = ""))
        
        embryo_overview <- data.frame(merged_CNVs_embryo_post)
        embryo_overview <- merge(embryo_overview, chr_lengths, all.x = T)
        embryo_overview$SV_size <- round(embryo_overview$width/embryo_overview$chr_size, 2)
        embryo_overview$CNV_ID <- paste(embryo_overview$seqnames, embryo_overview$start, embryo_overview$end, sep = "_")
        
        i <- 0
        
        print(paste("Counting CNVs for each cell in embryo: ",embryo ,sep = ""))
        for(sample in samplesheet$Cell[samplesheet$Embryo == embryo]){
          print(sample)
          sample_info <- samplesheet[which(samplesheet$Cell == sample),]
          
          if(sample_info$total.read.count > 10000){
            
            aneufinder_data_cell <- aneufinder_data[[names(aneufinder_data)[grep(paste(sample, ".bam", sep = ""), x = names(aneufinder_data))]]]$segments
            
            if(!is.null(aneufinder_data_cell)){
              
              cnvs_cell <- aneufinder_data_cell[which(aneufinder_data_cell$copy.number != sample_info$base_cn_state),]
              
              # Fragmented cells
              if(sample_info$base_cn_state == 0){
                cnvs_cell$copy.number <- 1
              }
              
              cnvs_cell <- Chain_CNVs(cnvs_cell)
              
              cnvs_cell$counts <- ifelse(cnvs_cell$copy.number > sample_info$base_cn_state, 1, -1)
              cnvs_cell$counts <- ifelse(cnvs_cell$copy.number > sample_info$base_cn_state+1, 2, cnvs_cell$counts)
              cnvs_cell$counts <- ifelse(cnvs_cell$copy.number > sample_info$base_cn_state+2, 3, cnvs_cell$counts)
              
              cnvs_cell$counts <- ifelse(cnvs_cell$copy.number < sample_info$base_cn_state-1, -2, cnvs_cell$counts)
              
              cnv_count_cell <- data.frame(Count_Overlapping_CNVs(x = merged_CNVs_embryo, y = cnvs_cell, minimal_overlap = 0.6, maximum_overlap = 0.9))
              cnv_count_cell$CNV_ID <- paste(cnv_count_cell$seqnames, cnv_count_cell$start, cnv_count_cell$end, sep = "_")
              
              embryo_overview <- merge(embryo_overview, cnv_count_cell[,c("CNV_ID", "counts")], all.x = T, by = "CNV_ID")
              names(embryo_overview)[ncol(embryo_overview)] <- sample
              embryo_overview <- embryo_overview[order(embryo_overview$seqnames, embryo_overview$start),]
              i <- i + 1
            }
          }
        }
      }
      
      if(nrow(embryo_overview)>0){
        embryo_overview[is.na(embryo_overview)] <- 0
        embryo_overview$Embryo_counts <- 0
        
        # Count how many each merged CNV occurs within the embryo
        for(i in 1:nrow(embryo_overview)){
          #print(embryo_overview[i,])
          embryo_overview$Embryo_counts[i] <- length(which(embryo_overview[i,-c(1:8)] != 0))
        }
        
        ## Remove the events that did not pass the filter settings
        # The coordinates of some events may have been changed due to chaining and merging etc
        unique_CNVS_filter_pass <- c()
        for(cell in names(embryo_overview)[-c(1:8, ncol(embryo_overview))]){
          #print(cell)
          unique_events_cell <- embryo_overview[which(embryo_overview[,cell] != 0 & embryo_overview$Embryo_counts == 1),]
          unique_events_cell_g <- GRanges(seqnames = unique_events_cell$seqnames, IRanges(start = unique_events_cell$start, end = unique_events_cell$end))
          #print(unique_events_cell)
          if(nrow(unique_events_cell) > 0 ){
            if(length(CNVs_filter_pass[[cell]]) > 0){
              #overlaps_unique_filter_pass <- findOverlaps(unique_events_cell_g, CNVs_filter_pass[[cell]])
              
              overlaps_unique_filter_pass <- Count_Overlapping_CNVs(unique_events_cell_g, CNVs_filter_pass[[cell]], minimal_overlap = 0.5, maximum_overlap = 0.9)
              filter_pass_IDs <- paste(seqnames(overlaps_unique_filter_pass), start(overlaps_unique_filter_pass), end(overlaps_unique_filter_pass), sep= "_")
              
              unique_CNVS_filter_pass <- c(unique_CNVS_filter_pass, unique(filter_pass_IDs))

            }
          }
        }
        if(length(unique_CNVS_filter_pass) > 0 ){
          embryo_overview <- embryo_overview[which(embryo_overview$CNV_ID %in% unique_CNVS_filter_pass | embryo_overview$Embryo_counts > 1),]
        } else {
          embryo_overview <- embryo_overview[which(embryo_overview$Embryo_counts > 1),]
        }
      }
      
      if(nrow(embryo_overview)>0){
        embryo_overview$Embryo_freq <- round( embryo_overview$Embryo_counts / (ncol(embryo_overview) -  9),2)
        
        cells <- ncol(embryo_overview) - 10
        embryo_overview$Cells_in_embryo <- cells
        
        for(i in 1:nrow(embryo_overview)){
          embryo_overview$Gains[i] <- length(which(embryo_overview[i,9:(8+cells)] > 0)) 
          embryo_overview$Losses[i] <- length(which(embryo_overview[i,9:(8+cells)] < 0)) 
        }
        
        # Reciprocal gains/losses should have a ~50:50 distribution (now: gains / losses > 30%)
        embryo_overview$Reciprocal <- ifelse(abs(embryo_overview$Gains / embryo_overview$Losses) > 0.1 &
                                           !is.infinite(abs(embryo_overview$Gains / embryo_overview$Losses)), "Reciprocal", "Non-reciprocal")
        
        # # The originating cell division is based on the embryo frequency, which should be higher than a dynamic threshold (depending on the number of sequenced cells in the embryo) 
        # thresholds_embryo <- threshold[[cells]]
        # 
        # embryo_overview$Division <- 4
        # 
        # if(cells > 3) {
        #   for(i in 1:nrow(embryo_overview)){
        #     
        #     embryo_overview$Division[i] <- ifelse(embryo_overview$Reciprocal[i] == "Non-reciprocal", 
        #                                       as.vector(thresholds_embryo$Division)[which(thresholds_embryo$Non_reciprocal_min <= embryo_overview$Embryo_freq[i] & 
        #                                                                                     thresholds_embryo$Non_reciprocal_max >= embryo_overview$Embryo_freq[i])],
        #                                       as.vector(thresholds_embryo$Division)[which(thresholds_embryo$Reciprocal_min <= embryo_overview$Embryo_freq[i] & 
        #                                                                                     thresholds_embryo$Reciprocal_max >= embryo_overview$Embryo_freq[i])])
        #   }
        # } else {
        #   embryo_overview$Division <- "ND"
        # }
        
        embryo_overview <- embryo_overview[!duplicated(embryo_overview$CNV_ID),]
        
        # Write the CNVs of the embryo:
        write.table(x = embryo_overview, file = paste(output_folder, "CNVs_", embryo, ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
        
        embryo_counts <- data.frame()
        
        # Count the number of CNVs per cell:
        for(cell in names(embryo_overview)[9:(9+unique(embryo_overview$Cells_in_embryo)-1)]){
          
          cnvs_cell <- embryo_overview[,c("CNV_ID","seqnames","start","end","width","strand","chr_size", "SV_size", "Embryo_counts","Reciprocal", names(embryo_overview)[names(embryo_overview) == cell])]
          cnvs_cell <- cnvs_cell[cnvs_cell[,cell] != 0,]
          
          if(nrow(cnvs_cell) > 0){
            cnvs_cell_g <- GRanges(seqnames = cnvs_cell$seqnames, IRanges(start = cnvs_cell$start, end = cnvs_cell$end), ID = cnvs_cell$CNV_ID)
            
            aneufinder_data_cell <- aneufinder_data[[names(aneufinder_data)[grep(paste(cell, ".bam", sep = ""), x = names(aneufinder_data))]]]$bins
            
            # determine which bins overlap with the filtered CNV calls:
            olap_CNVs_bins <- findOverlaps(cnvs_cell_g, aneufinder_data_cell)
            bins_per_CNVs <- cbind(as.data.frame(cnvs_cell[queryHits(olap_CNVs_bins),]), as.data.frame(aneufinder_data_cell[subjectHits(olap_CNVs_bins),]))
            
            # Calculate the median counts per CNV call
            median_counts <- aggregate(bins_per_CNVs$counts ~ bins_per_CNVs$CNV_ID, FUN = "median")
            names(median_counts) <- c("CNV_ID", "Median_counts")
            
            # Calculate the median counts for the bins with the base cn state (without filtered CNV calls)
            median_count_base_state <- median(aneufinder_data_cell[-subjectHits(olap_CNVs_bins),]$counts)
            
            cnvs_cell <- merge(cnvs_cell, median_counts, by = "CNV_ID", all.x = T)
            cnvs_cell$median_count_base_state <- median_count_base_state
            cnvs_cell$Relative_CNV_Change <- cnvs_cell$Median_counts/cnvs_cell$median_count_base_state
            
            dir.create( paste(input_folder, "Data/CNVs_cells/",resolution, "/", run, sep = ""), showWarnings = F, recursive = T)
            write.table(x = cnvs_cell, file = paste(input_folder, "Data/CNVs_cells/",resolution, "/", run, "/CNVs_", cell, ".txt", sep = ""), sep = "\t", quote = F, row.names = F)
          }
          
          if(nrow(cnvs_cell) > 0){
            cell_counts <- data.frame(Type = c("CNV", "CNV", 
                                               "CNV","CNV","CNV", 
                                               "Aneuploidy", "Aneuploidy",
                                               "Aneuploidy", "Aneuploidy", "Aneuploidy", "Filter", "Filter", "Filter"), 
                                      Ploidy = c("Loss", "Gain", 
                                                 "Recurrent", "Reciprocal", "Non-reciprocal",
                                                 "Loss", "Gain",
                                                 "Recurrent", "Reciprocal", "Non-reciprocal",  "Loss", "Gain", "Nullisomies"), 
                                      Counts = c(nrow(cnvs_cell[which(cnvs_cell$SV_size < 0.95 & cnvs_cell$width > 10e6 & cnvs_cell[,cell] < 0),]),
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size < 0.95 & cnvs_cell$width > 10e6 & cnvs_cell[,cell] > 0),]), 
                                                 
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size < 0.95 & cnvs_cell$width > 10e6 & 
                                                                        cnvs_cell$Embryo_counts > 1 & cnvs_cell$Reciprocal != "Reciprocal"),]),
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size < 0.95 & cnvs_cell$width > 10e6  & 
                                                                        cnvs_cell$Embryo_counts > 1 & cnvs_cell$Reciprocal == "Reciprocal"),]),
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size < 0.95 & cnvs_cell$width > 10e6  & 
                                                                        cnvs_cell$Embryo_counts == 1),]),
                                                 
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size > 0.95 & cnvs_cell$seqnames != "X" & cnvs_cell[,cell] < 0),]), 
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size > 0.95 & cnvs_cell$seqnames != "X" & cnvs_cell[,cell] > 0),]),
                                                 
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size > 0.95 & cnvs_cell$seqnames != "X" & 
                                                                        cnvs_cell$Embryo_counts > 1 & cnvs_cell$Reciprocal != "Reciprocal"),]),
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size > 0.95 & cnvs_cell$seqnames != "X" & cnvs_cell$Embryo_counts > 1 & cnvs_cell$Reciprocal == "Reciprocal"),]),
                                                 nrow(cnvs_cell[which(cnvs_cell$SV_size > 0.95 & cnvs_cell$seqnames != "X" & cnvs_cell$Embryo_counts == 1),]),
                                                 ifelse(is.null(filtered_events[[cell]]$Filtered_losses), 0, filtered_events[[cell]]$Filtered_losses),
                                                 ifelse(is.null(filtered_events[[cell]]$Filtered_gains), 0, filtered_events[[cell]]$Filtered_gains),
                                                 ifelse(is.null(filtered_events[[cell]]$Filter_nullisomies), 0, filtered_events[[cell]]$Filter_nullisomies)))
                                    
          } else {
            cell_counts <- data.frame(Type = c("CNV", "CNV", 
                                               "CNV","CNV","CNV", 
                                               "Aneuploidy", "Aneuploidy",
                                               "Aneuploidy", "Aneuploidy", "Aneuploidy", "Filter", "Filter", "Filter"), 
                                      Ploidy = c("Loss", "Gain", 
                                                 "Recurrent", "Reciprocal", "Non-reciprocal",
                                                 "Loss", "Gain",
                                                 "Recurrent", "Reciprocal", "Non-reciprocal",  "Loss", "Gain", "Nullisomies"), 
                                      Counts = c(0,0, 0,0,0,0, 0,0,0,0,
                                                 ifelse(is.null(filtered_events[[cell]]$Filtered_losses), 0, filtered_events[[cell]]$Filtered_losses),
                                                 ifelse(is.null(filtered_events[[cell]]$Filtered_gains), 0, filtered_events[[cell]]$Filtered_gains),
                                                 ifelse(is.null(filtered_events[[cell]]$Filter_nullisomies), 0, filtered_events[[cell]]$Filter_nullisomies)))
          }
            
          if(ncol(embryo_counts) > 0){
            embryo_counts <- cbind(embryo_counts, cell_counts[,ncol(cell_counts)])
          } else {
            embryo_counts <- cell_counts
          }
          names(embryo_counts)[ncol(embryo_counts)] <- cell
        }
        
        } else {
          print(paste("## No CNVs for embryo: ", embryo, sep = ""))
          embryo_counts <- data.frame()
          for(sample in samplesheet$Cell[samplesheet$Embryo == embryo]){
            cell_counts <- data.frame(Type = c("CNV", "CNV", 
                                               "CNV","CNV","CNV", 
                                               "Aneuploidy", "Aneuploidy",
                                               "Aneuploidy", "Aneuploidy", "Aneuploidy", "Filter", "Filter", "Filter"), 
                                      Ploidy = c("Loss", "Gain", 
                                                 "Recurrent", "Reciprocal", "Non-reciprocal",
                                                 "Loss", "Gain",
                                                 "Recurrent", "Reciprocal", "Non-reciprocal",  "Loss", "Gain", "Nullisomies"), 
                                      Counts = c(0,0, 0,0,0,0, 0,0,0,0,
                                                 ifelse(is.null(filtered_events[[sample]]$Filtered_losses), 0, filtered_events[[sample]]$Filtered_losses),
                                                 ifelse(is.null(filtered_events[[sample]]$Filtered_gains), 0, filtered_events[[sample]]$Filtered_gains),
                                                 ifelse(is.null(filtered_events[[sample]]$Filter_nullisomies), 0, filtered_events[[sample]]$Filter_nullisomies)))
            if(ncol(embryo_counts) > 0){
              embryo_counts <- cbind(embryo_counts, cell_counts[,ncol(cell_counts)])
            } else {
              embryo_counts <- cell_counts
            }
            names(embryo_counts)[ncol(embryo_counts)] <- sample
            
          }
          
        }
      write.table(x = embryo_counts, file = paste(CNV_Counts_Folder, "CNV_Counts_", embryo, ".txt", sep = ""), sep = "\t", quote = F, row.names = F)
    }
}
