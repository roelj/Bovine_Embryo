# This script is used to plot karyograms per embryo

library(AneuFinder)
library(ggplot2)
library(reshape2)

options(scipen = 999)

args <- commandArgs(trailingOnly=TRUE)

resolution <- args[1]
input_folder <- args[2]
plot_output_folder <- args[3]
runs <- args[4]

if(is.null(resolution)){
  resolution <- 1000000
} else {
  resolution <- as.numeric(resolution)
}

QC_folder <- paste(input_folder, "Data/QC/", sep = "")
samplesheet <-  read.delim(paste(QC_folder,"Overview_cells_",resolution,".txt", sep = ""), stringsAsFactors = F)
if(runs != ""){
  samplesheet <- samplesheet[which(samplesheet$Run %in% runs),]
}

dir.create(plot_output_folder)
aneufinder_folder <-  paste(input_folder,"Data/Aneufinder/", sep = "")
cnv_folder <-  paste(input_folder, "Data/CNV/", resolution, "/", sep = "")
cn_folder <-  paste(input_folder, "Data/CN_States/", "/", sep = "")
SNP_density_folder <- paste(input_folder, "Data/SNVs/SNP_Counts/Bins_10Mb/", sep = "")
SNP_FW_folder <- paste(input_folder, "Data/SNVs/SNP_Counts/Bins_10Mb_FW/", sep = "")
SNP_RV_folder <- paste(input_folder, "Data/SNVs/SNP_Counts/Bins_10Mb_RV/", sep = "")

## This function plots data from a single chromosome showing the read distribution and the copy number state
# The height of the karyogram of one cell is determined by max_height
plot_single_chromosome <- function(chromosome, 
                                   method,
                                   cell_data = data, 
                                   cnvs = NULL,
                                   breakpoints_chr = "",
                                   max_score, 
                                   start = 0,
                                   total_length, width = 1000, 
                                   y0 = 0,
                                   max_height = 2,
                                   states = 3, 
                                   background = TRUE,
                                   SNPs_density = "",
                                   empty = FALSE, 
                                   SNPs_FW = "",
                                   SNPs_RV = "",
                                   show_phasing = TRUE, base_cn_state = 2){
  
  chromosome_data <- as.data.frame(cell_data$bins[seqnames(cell_data$bins) == chromosome])
  #print(head(chromosome_data))
  # This is the height of the boxes showing the CN state and SNP densities
  box_height <- ifelse(show_phasing == TRUE, 0.25, 0.3)
  gap_between_boxes <-ifelse(show_phasing == TRUE, 0.075, 0)
  

  number_boxes <- ifelse(show_phasing == TRUE & method == "Strand-seq", 3, 1)
  number_boxes <- ifelse(show_phasing == TRUE & method != "Strand-seq", 2, number_boxes)
  
  if(empty == FALSE){
    if(method == "WGS"){
      #print(head(chromosome_data))
      chromosome_data$color <- ifelse(chromosome_data$copy.number == 2, "chartreuse3", "gray")
      chromosome_data$color <- ifelse(chromosome_data$copy.number == 1, "#E41A1C", chromosome_data$color)
      chromosome_data$color <- ifelse(chromosome_data$copy.number == 3, "#2A6AFF", chromosome_data$color)
      chromosome_data$color <- ifelse(chromosome_data$copy.number == 4, "#1034A6", chromosome_data$color)
      chromosome_data$color <- ifelse(chromosome_data$copy.number > 4,  "#111E6C", chromosome_data$color)
      
    } else {
      chromosome_data$mcolor <- ifelse(chromosome_data$mcopy.number == 2, "chartreuse3", "gray")
      chromosome_data$mcolor <- ifelse(chromosome_data$mcopy.number == 1, "#E41A1C", chromosome_data$mcolor)
      chromosome_data$mcolor <- ifelse(chromosome_data$mcopy.number == 3, "#2A6AFF", chromosome_data$mcolor)
      chromosome_data$mcolor <- ifelse(chromosome_data$mcopy.number == 4, "#1034A6", chromosome_data$mcolor)
      chromosome_data$mcolor <- ifelse(chromosome_data$mcopy.number > 4, "#111E6C", chromosome_data$mcolor)
      
      chromosome_data$pcolor <- ifelse(chromosome_data$pcopy.number == 2, "chartreuse3", "gray")
      chromosome_data$pcolor <- ifelse(chromosome_data$pcopy.number == 1, "#E41A1C", chromosome_data$pcolor)
      chromosome_data$pcolor <- ifelse(chromosome_data$pcopy.number == 3, "#2A6AFF", chromosome_data$pcolor)
      chromosome_data$pcolor <- ifelse(chromosome_data$pcopy.number == 4, "#1034A6", chromosome_data$pcolor)
      chromosome_data$pcolor <- ifelse(chromosome_data$pcopy.number > 4, "#111E6C", chromosome_data$pcolor)
    }

  }
  
  # This is for the background color of the chromosome 
  if(background == TRUE){
    rect(ybottom = y0, ytop = y0+max_height, 
         xleft = start + min(chromosome_data$start) / total_length * width, 
         xright = start + max(chromosome_data$end) / total_length * width, 
         border = NA, col = rgb(211,211,211, 100, maxColorValue = 255))
  }
  
  # This is to plot the coverage in bins:
  if(empty == FALSE){
    if(method == "WGS"){
      #states <- 2
      rect(ybottom = y0+box_height/2, 
           ytop = ifelse((chromosome_data$counts / max_score) * max_height < max_height, 
                         y0+box_height/2 + (chromosome_data$counts / max_score)  * max_height, 
                         y0+box_height/2 + max_height),
           xleft = start + chromosome_data$start / total_length * width, 
           xright = start + chromosome_data$end / total_length * width,
           border = NA, col = chromosome_data$color)
    } else if (method == "Strand-seq"){
      
      states <- 1.5
      # The mcounts will be plotted below the chromosome bar, the pcounts above the bar:
      chromosome_data$mcounts <- chromosome_data$mcounts*-1

      # For strand-seq data the y0 is in the center of the karyogram
      y0 <- y0+max_height/2
      
      y0_coverage <- ifelse(show_phasing == TRUE, 
                            box_height/2 + gap_between_boxes + box_height,
                            box_height/2)
      
      
      #y0+box_height/2 + gap_between_boxes, ytop = y0+box_height/2 + gap_between_boxes + box_height
      
      # Plot the pcounts (FW strand)
      rect(ybottom = y0+y0_coverage, 
           ytop = ifelse((abs(chromosome_data$pcounts / max_score) / states  * max_height) < max_height/2, 
                         y0+y0_coverage + (chromosome_data$pcounts / max_score) / states * max_height, 
                         y0+y0_coverage + max_height/2),
           xleft = start + chromosome_data$start / total_length * width, 
           xright = start + chromosome_data$end / total_length * width,
           border = NA, col = chromosome_data$pcolor)
      # Plot the mcounts (RV strand)
      rect(ybottom = y0+(y0_coverage*-1), 
           ytop = ifelse((abs(chromosome_data$mcounts / max_score) / states * max_height) < max_height/2, 
                         y0+(y0_coverage*-1) + (chromosome_data$mcounts / max_score) / states * max_height, 
                         y0+(y0_coverage*-1) + max_height/-2),
           xleft = start + chromosome_data$start / total_length * width, 
           xright = start + chromosome_data$end / total_length * width,
           border = NA, col = chromosome_data$mcolor)
    }


  # Select the filtered CNVs for the cell. If they are not available the "raw" CN states determined by Aneufinder will be used (mostly for cells < 100k reads). 
  if(!is.null(cnvs) == TRUE){

    base_cn_state <- as.numeric(base_cn_state)
    cnvs_cell <-  cnvs[cnvs$seqnames == chromosome,]
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number == 2, rgb(102,205,0, 255, maxColorValue = 255), "gray")
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number == 1, rgb(228, 26, 28, 255, maxColorValue = 255), cnvs_cell$color)
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number == 3, rgb(55, 126, 184, 255, maxColorValue = 255), cnvs_cell$color)
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number == 4, rgb(16, 52, 166, 255, maxColorValue = 255), cnvs_cell$color)
    
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number > 4, rgb(17,30,108, 255, maxColorValue = 255), cnvs_cell$color)
    
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number == 2 & cnvs_cell$Copy_number == base_cn_state, rgb(102,205,0, 120, maxColorValue = 255), cnvs_cell$color)
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number == 1 & cnvs_cell$Copy_number == base_cn_state, rgb(228, 26, 28, 120, maxColorValue = 255), cnvs_cell$color)
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number == 3 & cnvs_cell$Copy_number == base_cn_state, rgb(55, 126, 184, 120, maxColorValue = 255), cnvs_cell$color)
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number == 4 &cnvs_cell$Copy_number == base_cn_state,  rgb(16, 52, 166, 120, maxColorValue = 255), cnvs_cell$color)
    cnvs_cell$color <- ifelse(cnvs_cell$Copy_number > 4 & cnvs_cell$Copy_number == base_cn_state, rgb(17,30,108, 120, maxColorValue = 255), cnvs_cell$color)
    #print(cnvs_cell)
    
  }
  
  # Plot a box below each karyogram containing the base copy number state
  rect(ybottom = ifelse(method == "WGS", y0-box_height/2, y0-box_height/2), ytop = ifelse(method == "WGS", y0+box_height/2, y0+box_height/2), 
       xleft = start + min(chromosome_data$start) / total_length * width,
       xright = start + max(chromosome_data$end) / total_length * width,
       border = "white", col = "white", lwd = 1)
  
  # Plot the CNVs on top of this box with colors indicating the CN
  if(nrow(cnvs_cell) > 0){
    rect(ybottom = ifelse(method == "WGS", y0-box_height/2, y0-box_height/2), ytop = ifelse(method == "WGS",  y0+box_height/2, y0+box_height/2), 
         xleft = start + cnvs_cell$start / total_length * width,
         xright = start + cnvs_cell$end / total_length * width,
         border = "black", col = cnvs_cell$color)
  }


  # Plot the maternal snp density
  if(method == "WGS"){
    if(show_phasing == TRUE){
      if(length(SNPs_density) > 1){
        #print(SNPs_density)
        SNPs_density <- SNPs_density[SNPs_density$MAT_Counts + SNPs_density$PAT_MAT_Counts > 10,]
        
        if(nrow(SNPs_density) > 0){
          # Plot a box for the PAT/MAT SNP density
          rect(ybottom = y0-gap_between_boxes-box_height/2-box_height, 
               ytop = y0-gap_between_boxes-box_height/2, 
               xleft = start + min(chromosome_data$start) / total_length * width,
               xright = start + max(chromosome_data$end) / total_length * width,
               border = "black", col = "white", lwd = 1)
          
          SNPs_density$MAT_Freq[is.na(SNPs_density$MAT_Freq)] <- 1
          SNPs_density$color <- "white"
          SNPs_density$alpha <- 1
          SNPs_density$alpha[SNPs_density$MAT_Freq < 0.25] <- 0.3 + (0.15-SNPs_density$MAT_Freq[SNPs_density$MAT_Freq < 0.25])/0.2
          SNPs_density$alpha[SNPs_density$MAT_Freq > 0.5] <- 0.3 + (SNPs_density$MAT_Freq[SNPs_density$MAT_Freq > 0.5]-0.5)/0.5
          SNPs_density$alpha[SNPs_density$alpha > 1] <- 1
          SNPs_density$alpha[SNPs_density$alpha < 0] <- 0
          
          SNPs_density$color[SNPs_density$MAT_Freq < 0.25] <- rgb(255,128,0, SNPs_density$alpha[SNPs_density$MAT_Freq < 0.25]*255, maxColorValue = 255)
          SNPs_density$color[SNPs_density$MAT_Freq > 0.5] <- rgb(127,56,236,  SNPs_density$alpha[SNPs_density$MAT_Freq > 0.5]*255, maxColorValue = 255)
          
          #print( SNPs_density[SNPs_density$MAT_Freq > 0.4,])
          
          #SNPs_density$color <- ifelse(SNPs_density$MAT_Freq > 0.45 & SNPs_density$MAT_Freq < 1, "darkred", SNPs_density$color)
          rect(ybottom = y0-gap_between_boxes-box_height/2-box_height, 
               ytop = y0-gap_between_boxes-box_height/2, 
               xleft = start + SNPs_density$start / total_length * width,
               xright = start + SNPs_density$end / total_length * width,
               border = NA, col = SNPs_density$color)
          
          rect(ybottom = y0-gap_between_boxes-box_height/2-box_height, 
               ytop = y0-gap_between_boxes-box_height/2, 
               xleft = start + min(chromosome_data$start) / total_length * width,
               xright = start + max(chromosome_data$end) / total_length * width,
               border = "black", col = NA, lwd = 1)
        }
 
  } else {
    # Plot a box for the PAT/MAT SNP density
    rect(ybottom = y0-gap_between_boxes-box_height/2-box_height, 
         ytop = y0-gap_between_boxes-box_height/2,  
         xleft = start + min(chromosome_data$start) / total_length * width,
         xright = start + max(chromosome_data$end) / total_length * width,
         border = "black", col = "gray", lwd = 1)
  }
    }
    
    
    } else if (method == "Strand-seq"){
      if(show_phasing == TRUE){
        if(length(SNPs_FW) > 1 | length(SNPs_RV) > 1){
          
          SNP_density_FW <- SNPs_FW[which(SNPs_FW$MAT_Counts + SNPs_FW$PAT_MAT_Counts > 10),]
          SNP_density_RV <- SNPs_RV[which(SNPs_RV$MAT_Counts + SNPs_RV$PAT_MAT_Counts > 10),]
          
          if(nrow(SNP_density_FW) > 0){
            
            # background for SNP density track FW strand
            rect(ybottom =y0+box_height/2 + gap_between_boxes, ytop = y0+box_height/2 + gap_between_boxes + box_height, 
                 xleft = start + min(chromosome_data$start) / total_length * width,
                 xright = start + max(chromosome_data$end) / total_length * width,
                 border = "black", col = "white", lwd = 1)
            
            SNP_density_FW$MAT_Freq[is.na(SNP_density_FW$MAT_Freq)] <- 1
            SNP_density_FW$color <- "white"
            SNP_density_FW$alpha <- 1
            SNP_density_FW$alpha[SNP_density_FW$MAT_Freq < 0.25] <- 0.3 + (0.15-SNP_density_FW$MAT_Freq[SNP_density_FW$MAT_Freq < 0.25])/0.2
            SNP_density_FW$alpha[SNP_density_FW$MAT_Freq > 0.5] <- 0.3 + (SNP_density_FW$MAT_Freq[SNP_density_FW$MAT_Freq > 0.5]-0.5)/0.5
            SNP_density_FW$alpha[SNP_density_FW$alpha > 1] <- 1
            SNP_density_FW$alpha[SNP_density_FW$alpha < 0] <- 0
            SNP_density_FW$color[SNP_density_FW$MAT_Freq < 0.25] <- rgb(255,128,0, SNP_density_FW$alpha[SNP_density_FW$MAT_Freq < 0.25]*255, maxColorValue = 255)
            SNP_density_FW$color[SNP_density_FW$MAT_Freq > 0.5] <- rgb(127,56,236,  SNP_density_FW$alpha[SNP_density_FW$MAT_Freq > 0.5]*255, maxColorValue = 255)
            
            # Plot the mat snp density
            rect(ybottom = y0+box_height/2 + gap_between_boxes, ytop = y0+box_height/2 + gap_between_boxes + box_height, 
                 xleft = start + SNP_density_FW$start / total_length * width,
                 xright = start + SNP_density_FW$end / total_length * width,
                 border = NA, col = SNP_density_FW$color)
            rect(ybottom = y0+box_height/2 + gap_between_boxes, ytop = y0+box_height/2 + gap_between_boxes + box_height, 
                 xleft = start + min(chromosome_data$start) / total_length * width,
                 xright = start + max(chromosome_data$end) / total_length * width,
                 border = "black", col = NA, lwd = 1)
          }
          
          
          if(nrow(SNP_density_RV) > 0){
            # plot background
            # y0+box_height/2 + gap_between_boxes, ytop = y0+box_height/2 + gap_between_boxes + box_height
            # y0_coverage <- box_height/2 + gap_between_boxes + box_heigh (start coverage) 
            # y0+(y0_coverage*-1)
            rect(ybottom = y0-(box_height/2+gap_between_boxes+box_height), ytop = y0-(box_height/2+gap_between_boxes), 
                 xleft = start + min(chromosome_data$start) / total_length * width,
                 xright = start + max(chromosome_data$end) / total_length * width,
                 border = "black", col = "white", lwd = 1)
            
            SNP_density_RV$MAT_Freq[is.na(SNP_density_RV$MAT_Freq)] <- 1
            SNP_density_RV$color <- "white"
            SNP_density_RV$alpha <- 1
            SNP_density_RV$alpha[SNP_density_RV$MAT_Freq < 0.25] <- 0.3 + (0.15-SNP_density_RV$MAT_Freq[SNP_density_RV$MAT_Freq < 0.25])/0.2
            SNP_density_RV$alpha[SNP_density_RV$MAT_Freq > 0.5] <- 0.3 + (SNP_density_RV$MAT_Freq[SNP_density_RV$MAT_Freq > 0.5]-0.5)/0.5
            SNP_density_RV$alpha[SNP_density_RV$alpha > 1] <- 1
            SNP_density_RV$alpha[SNP_density_RV$alpha < 0] <- 0
            SNP_density_RV$color[SNP_density_RV$MAT_Freq < 0.25] <- rgb(255,128,0, SNP_density_RV$alpha[SNP_density_RV$MAT_Freq < 0.25]*255, maxColorValue = 255)
            SNP_density_RV$color[SNP_density_RV$MAT_Freq > 0.5] <- rgb(127,56,236,  SNP_density_RV$alpha[SNP_density_RV$MAT_Freq > 0.5]*255, maxColorValue = 255)
            
            rect(ybottom = y0-(box_height/2+gap_between_boxes+box_height), ytop = y0-(box_height/2+gap_between_boxes), 
                 xleft = start + SNP_density_RV$start / total_length * width,
                 xright = start + SNP_density_RV$end / total_length * width,
                 border = NA, col = SNP_density_RV$color)
            # plot the border around mat SNP density track
            rect(ybottom = y0-(box_height/2+gap_between_boxes+box_height), ytop = y0-(box_height/2+gap_between_boxes), 
                 xleft = start + min(chromosome_data$start) / total_length * width,
                 xright = start + max(chromosome_data$end) / total_length * width,
                 border = "black", col = NA, lwd = 1)
          }
        }
      }
    }
  } else {
    # plot a black box if the cells has less reads than the threshold of 10k reads
    rect(ybottom = y0-0.1, ytop = y0, 
         xleft = start + min(chromosome_data$start) / total_length * width,
         xright = start + max(chromosome_data$end) / total_length * width,
         border = "black", col = "black", lwd = 1)
  }
}



## This function plots multiple single chromosomes in one row
# This funtions mainly determines the widths of each chromosome
# It also determines the max_score value (median count of bins with copy number state 2 times 1.5)
plot_multiple_chromosomes <- function(y0, input_data, 
                                      cnvs,
                                      gap_between_chr = 5, 
                                      plot_width = 1000, 
                                      chr_label = FALSE, 
                                      chromosomes = c(1:29, "X"), 
                                      cell_height = 2,
                                      empty = FALSE, method,
                                      breakpoints = "",
                                      SNPs_cell = "", 
                                      SNPs_FW = "",
                                      SNPs_RV = "", 
                                      show_phasing = TRUE,
                                      number_cells = 8,
                                      Ploidy = "Diploid",
                                      base_cn_state = 2){
  
  # The plot has a fixed width (1000) on the X-axis. This width is divided over the chromosomes and the gaps between the chromosomes.
  # The width factor is the total width available for the chromosomes (total plot width - width of all the gaps)
  width_factor <- plot_width - gap_between_chr * length(chromosomes)

  cn_states_cell <- input_data$bins
  
  overview <- data.frame()
  
  # Calculate the total width (in basepairs) for all chromosomes that are plotted together:
  total_width_bp <- 0
  for(chr in unique(seqnames(input_data$bins))){
    total_width_bp <- total_width_bp+max(end(input_data$bins[seqnames(input_data$bins) == chr]))
  }
  
  # Determine the plot start and end positions of each chromosome (based on the basepair width of the chr and the plot width)
  for(chr in unique(seqnames(input_data$bins))){
    #print(chr)
    
    if(chr != "1"){
      chr_overview <- data.frame(chr = chr, start = overview[nrow(overview), "end"] + gap_between_chr, 
                                 end = overview[nrow(overview), "end"] + gap_between_chr + max(end(input_data$bins[seqnames(input_data$bins) == chr])) /  total_width_bp * width_factor)
      overview <- rbind(overview, chr_overview)
     } else {
      chr_overview <- data.frame(chr = chr, start = 0, end = (0+max(end(input_data$bins[seqnames(input_data$bins) == chr]))) / total_width_bp * width_factor)
      overview <- rbind(overview, chr_overview)
    }
  }
  
  # Alter copy number state for haploids to 1
  if(Ploidy == "Haploid"){
    #input_data$segments$copy.number <- ceiling(input_data$segments$copy.number / 2)
    if(length(which(input_data$bins$copy.number == 1)) < 100){
      input_data$bins$copy.number <- ceiling(input_data$bins$copy.number / 2)
    }
  }
  
  
  # Plot all the chromosomes one by one
  for(chr in unique(seqnames(input_data$bins))){
    #print(chr)
    if(chr != "X"){
      # Uneven chromosomes get a background (X is not even or uneven and therefore excluded)
      background <- ifelse(as.numeric(chr) %% 2 == 0, FALSE, TRUE )
    }
    
    if(empty == FALSE){
      
      max_score <- as.numeric(quantile(cn_states_cell$counts, 0.95))*2
      max_state <- max(cn_states_cell$copy.number[cn_states_cell$counts < max_score])

      if(length(SNPs_cell) > 1){
        SNPs_chr <- SNPs_cell[which(SNPs_cell$seqnames == chr),]
      } else {
        SNPs_chr <- ""
      }
      
      if(length(SNPs_FW) > 1){
        SNPs_FW_chr <- SNPs_FW[which(SNPs_FW$seqnames == chr),]
      } else {
        SNPs_FW_chr <- ""
      }
      if(length(SNPs_RV) > 1){
        SNPs_RV_chr <- SNPs_RV[which(SNPs_RV$seqnames == chr),]
      } else {
        SNPs_RV_chr <- ""
      }
      plot_single_chromosome(chromosome = chr, 
                             method = method,
                             start = overview$start[overview$chr == chr], 
                             total_length = total_width_bp, 
                             width = width_factor, y0 = y0, 
                             cell_data = input_data, 
                             cnvs = cnvs,
                             max_score = max_score,
                             background = background, max_height = cell_height, 
                             SNPs_density = SNPs_chr,
                             SNPs_FW = SNPs_FW_chr,
                             SNPs_RV = SNPs_RV_chr,
                             show_phasing = show_phasing, base_cn_state = base_cn_state)
      
    } else { 
      # Plot empty cell (containing only a black bar per chr and the background)
      plot_single_chromosome(chromosome = chr, 
                             method = method,
                             start = overview$start[overview$chr == chr], 
                             total_length = total_width_bp, 
                             width = width_factor, y0 = y0, 
                             cell_data = input_data, 
                             cnvs = cnvs,
                             max_score = median(cn_states_cell$counts[cn_states_cell$copy.number == 2])*1.5,
                             background = background, max_height = cell_height,
                             empty = TRUE,
                             show_phasing = show_phasing)
    }
  }
  
  if(chr_label == TRUE){
    if(number_cells > 2){
      text(y = y0-0.6, x = (overview$start+overview$end)/2, labels = overview$chr, cex = 0.8)
      
    } else {
      text(y = y0+3.75, x = (overview$start+overview$end)/2, labels = overview$chr, cex = 0.8, font = 2)
      
    }
  }
  
}


# This function generates the entire plot
# It determines the height of each cell based on the total height of the plot and the number of cells that have to be plotted.
plot_embryo <- function(embryo, 
                        run, 
                        aneufinder_data, 
                        cnvs, 
                        samplesheet = QC_data, 
                        gap_between_cells = 0.6, 
                        method, 
                        ymax = 20, 
                        breakpoints = "",
                        show_phasing = TRUE, cn_folder){
  
  # Select the embryo metadata 
  embryo_info <- samplesheet[which(samplesheet$Embryo == embryo & samplesheet$Run == run),]
  
  if(!is.null(embryo_info$Position)){
    # Reorder the metadata based on the position determined by cell clustering:
    embryo_info <- embryo_info[order(embryo_info$Position, decreasing = T),]
  }
  
  cells_to_plot <- nrow(embryo_info[which(embryo_info$total.read.count > 10000),])
  
  # Determine the height of each cell (on y-axis). Normal cells have 2x this height and empty cells without reads will have 1x this value
  # 
  height_per_cell <- ymax / (cells_to_plot*2 + nrow(embryo_info)-cells_to_plot + gap_between_cells*(nrow(embryo_info)-1))
  
  # The y_cumulative will determine where the next cell will be plotted. The first cell will be plotted at the ymax - the height of that cell
  y_cumulative <- ymax
  
  for(i in 0:(nrow(embryo_info)-1)){
    cell <- embryo_info$Cell[i+1]
    print(cell)
    
    cnvs <- read.delim(paste(cn_folder, "CN_States_", cell, ".txt", sep = ""), stringsAsFactors = F)
    
    chr_label <- ifelse(i == (nrow(embryo_info)-1), TRUE, FALSE)
    
    # Only plot cells with > 10000 reads, else plot an "empty" cell (just a black bar per chr)
    if(embryo_info[embryo_info$Cell == cell, "total.read.count"] > 10000){
      
      if(embryo_info[embryo_info$Cell == cell, "total.read.count"] > 10000){
        if(method == "WGS"){
          SNPs_cell <- read.delim(paste(SNP_density_folder, run, "/", cell, ".txt", sep = ""), stringsAsFactors = F)
          SNPs_FW <- ""
          SNPs_RV <- ""
        } else if (method == "Strand-seq"){
          SNPs_FW <- read.delim(paste(SNP_FW_folder, run, "/", cell, ".txt", sep = ""), stringsAsFactors = F)
          SNPs_RV <- read.delim(paste(SNP_RV_folder, run, "/", cell, ".txt", sep = ""), stringsAsFactors = F)
          SNPs_cell <- ""
        }
      } else {
        SNPs_cell <- ""
        SNPs_FW <- ""
        SNPs_RV <- ""
      }
      
      # y-cell is the bottom y-coordinate of the cell
      y_cell <- y_cumulative - height_per_cell*2
      
      Cell_ploidy <- embryo_info[embryo_info$Cell == cell, "Ploidy"]
      base_cn_state <- embryo_info[embryo_info$Cell == cell, "base_cn_state"]
      plot_multiple_chromosomes(y0 = y_cell, 
                                input_data = aneufinder_data[[grep(pattern = paste(cell, ".bam", sep = ""), x = names(aneufinder_data))]], 
                                chr_label = chr_label,
                                cell_height = height_per_cell*2,
                                cnvs = cnvs,
                                method = method,
                                SNPs_cell = SNPs_cell,
                                SNPs_FW = SNPs_FW, 
                                SNPs_RV = SNPs_RV, show_phasing = show_phasing,
                                number_cells = cells_to_plot, Ploidy = Cell_ploidy, base_cn_state = base_cn_state)
      
      # Plot the labels at the Y-axis
      # Cell ID
      
      
      Cell_ID <- gsub(pattern = paste(embryo_info$Run[embryo_info$Cell == cell], "_", sep = ""), replacement = "", x = cell)
      mtext(text = paste(embryo_info$Cell[which(embryo_info$Cell == cell)], Cell_ID, sep = " / "), side = 2, line = -2, at = y_cell+height_per_cell*1.5, cex = 0.8, las = 1,
            col = ifelse(embryo_info$Exclude[which(embryo_info$Cell == cell)] != "No", "darkred", "black")) 
      # Number of reads divided by 1 million
      mtext(text = paste(as.character(round(embryo_info[which(embryo_info$Cell == cell),"total.read.count"]/1e6,3)),"M reads", sep = ""), 
            side = 2, line = -2, at = y_cell+height_per_cell, cex = 0.8, las = 1, col = ifelse(embryo_info$Exclude[which(embryo_info$Cell == cell)] != "No", "darkred", "black"))
      
      mtext(text = as.character(embryo_info[which(embryo_info$Cell == cell),"Sex"]), side = 2, line = -2, at = y_cell+height_per_cell*0.5, cex = 0.8, las = 1)
      
      # conclusion
      mtext(text = as.character(embryo_info[which(embryo_info$Cell == cell),"Conclusion"]), side = 2, line = -2, at = y_cell, cex = 0.8, las = 1)
      
      y_cumulative <- y_cell-gap_between_cells
    } else {
      y_cell <- y_cumulative - height_per_cell
      
      plot_multiple_chromosomes(y0 = y_cell, 
                                input_data = aneufinder_data[[grep(pattern = paste(cell, ".bam", sep = ""), x = names(aneufinder_data))]], 
                                chr_label = chr_label,
                                cell_height = height_per_cell,
                                cnvs = cnvs,
                                empty = TRUE,
                                method = method, show_phasing = show_phasing, number_cells = cells_to_plot)
      Cell_ID <- gsub(pattern = paste(embryo_info$Run[embryo_info$Cell == cell], "_", sep = ""), replacement = "", x = cell)
      mtext(text = paste(embryo_info$Cell[which(embryo_info$Cell == cell)], Cell_ID, sep = " / "), side = 2, line = -2, at = y_cell+height_per_cell/1.5, cex = 0.8, las = 1,
            col = ifelse(embryo_info$Exclude[which(embryo_info$Cell == cell)] != "No", "darkred", "black")) 
      mtext(text = paste(as.character(round(embryo_info[which(embryo_info$Cell == cell),"total.read.count"]/1e6,3)),"M reads", sep = ""), side = 2, line = -2, at = y_cell+height_per_cell/3.5, cex = 0.8, las = 1)
      
      y_cumulative <- y_cell-gap_between_cells
    }
  }
  
}

show_phasing = TRUE
#show_phasing = FALSE

for(run in unique(samplesheet$Run)){
  print("# Reading Aneufinder data")
  print(run)
  
  aneufinder_output <- ifelse(length(grep("Strand-seq", x = unique(samplesheet[samplesheet$Run == run, "Method"])) > 0),
                              paste(aneufinder_folder, run, "/MODELS-StrandSeq/method-edivisive/", sep = ""),
                              paste(aneufinder_folder, run, "/MODELS/method-edivisive/", sep = ""))
  aneufinder_files <- list.files(aneufinder_output, full.names = T)

  # Only read files for the right resolution (Aneufinder folders may contain multiple datasets for different binsizes)
  if(length(grep(pattern = paste("bam_binsize_", format(resolution, scientific = T), sep = ""), x = as.character(aneufinder_files), fixed = T)) > 0){
    aneufinder_data <- loadFromFiles(aneufinder_files[grep(pattern = paste("bam_binsize_", format(resolution, scientific = T), sep = ""), x = as.character(aneufinder_files), fixed = T)])
    for(embryo in unique(samplesheet$Embryo[samplesheet$Run == run ])){
      print(embryo)
      
      embryo_stats <- samplesheet[which(samplesheet$Embryo == embryo & startsWith(x =  samplesheet$Embryo, prefix = "E")),]
      cnv_file <- paste(cnv_folder, "CNVs_", embryo, ".txt", sep = "" )
      cn_folder_run <- paste(cn_folder, run, "/", sep = "")
      if(file.exists(cnv_file) & file.info(cnv_file)$size > 1){
        cnvs <- read.delim(cnv_file, stringsAsFactors = F)
      } else {
        cnvs <- ""
      }
      
      if(nrow(embryo_stats) > 0){
        dir.create(paste(plot_output_folder, run, sep =""))
        
        output_height <- ifelse(unique(embryo_stats$Stage == "8-cell"), 7, 2.5)
        
        pdf(file = paste(plot_output_folder, run, "/", embryo, ".pdf", sep = ""), width = 10, height = output_height, pointsize = 10, family =  "ArialMT")
        
        par(mar=c(0, 4, 0, 0))
        plot.new()  
        
        ymax <- ifelse(unique(embryo_stats$Stage == "8-cell"), 20, 8)
        #print(ymax)
        xlim <- c(0,1000)
        ylim <- c(0, ymax)
        
        plot.window(xlim=xlim, ylim=ylim)
        title(main = paste(unique(samplesheet$Embryo[samplesheet$Run == run & samplesheet$Embryo == embryo]),
                           " (",unique(samplesheet$Group[samplesheet$Run == run & samplesheet$Embryo == embryo]), " Gy, ", 
                           unique(samplesheet$Embryo_classification[samplesheet$Run == run & samplesheet$Embryo == embryo]),
                           ")", sep = ""), line = -1)
        plot_embryo(embryo = embryo, run = run, aneufinder_data = aneufinder_data, cnvs = cnvs, samplesheet = embryo_stats, 
                    method = ifelse(length(grep("Strand-seq", x= embryo_stats$Method )) > 0, "Strand-seq", "WGS"), 
                    ymax = ymax, breakpoints = breakpoints, show_phasing = show_phasing, cn_folder = cn_folder_run)
        dev.off() 
      } else {
        print("# Embryo data not found")
      }
    }
  } else {
    print(paste("# Aneufinder data run ", run, " for resolution ", resolution, " NOT found!", sep = ""))
  }
}






