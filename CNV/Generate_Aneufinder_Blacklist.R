# This script can be used to generate a blacklist for aneufinder based on genomic coverage calculated from a BAM file
# Samtools idxstats will be used to determine the chromosome ordering in the bam file
args <- commandArgs(trailingOnly=TRUE)
options(scipen = 999)

library(GenomicRanges)
library(ggplot2)
library(AneuFinder)
BAM_file <- args[1]
output_folder <- args[2] # counts per bin will be written to this folder
blacklist_file <- args[3] # blacklist will be written to this file (in the outputfolder)
bin_size <- args[4] # coverage will be calculated for this bin size (100000 bp recommended by Aneufinder manual)

print(paste("# Generating blacklist (binsize = ",bin_size,") based on ", BAM_file, sep =""))

if(dir.exists(output_folder) == FALSE){
  dir.create(output_folder, recursive = T)
}

if(file.exists(BAM_file) == FALSE){
  stop(paste("! BAM file not found: ", BAM_file, sep =""))
}


# Obtain chr lengths from BSgenome.Btaurus.UCSC.bosTau8
# Remove unwanted chromosomes
autosomes <- 1:29
chromosomes <- c(1:29, "M", "X")

genome_file <- paste(output_folder, blacklist_file, ".idxstats", sep = "")

# Obtain the ordering of the chromosomes in the bam file using samtools idxstats:
samtools_command <- paste("samtools idxstats ", BAM_file, " > ", genome_file)
print(samtools_command)
system(samtools_command)

# Read the idxstats file
chr_sizes <- read.delim(genome_file, header = F, stringsAsFactors = F)
chr_sizes_g <- GRanges(seqnames = chr_sizes[which(chr_sizes[,1] %in% chromosomes),1], IRanges(start = 1, end = chr_sizes[which(chr_sizes[,1] %in% chromosomes),2]))
seqlengths(chr_sizes_g) <- chr_sizes[which(chr_sizes[,1] %in% chromosomes),2]

# Generate genomic bins with length of "bin_size"
bins <- tileGenome(seqlengths = seqlengths(chr_sizes_g), 
                   tilewidth = as.numeric(bin_size), cut.last.tile.in.chrom=TRUE)
seqlevels(bins) <- gsub(seqlevels(bins), pattern = "chr", replacement = "")
seqlevels(bins)[seqlevels(bins) == "M"] <- "MT"

bins_output <- data.frame(bins, stringsAsFactors = F)[,1:3]

# write the bins
bins_file <- paste(output_folder, "Genomic_Bins_", bin_size, "bp.bed", sep = "")
write.table(bins_output, file = bins_file, sep = "\t", quote = F, col.names = F, row.names = F)

# Run bedtools to calculate the coverage per bin
bedtools_output <- paste(output_folder, "Genome_Coverage_Embryos.bed", sep = "")
bedtools_command <- paste("bedtools intersect -a ",bins_file," -b ",BAM_file," -sorted -c -g ",genome_file," > ", bedtools_output,sep = "")
print("# Running bedtools to calculate coverage, may take a while")
print(bedtools_command)
system(bedtools_command)

print(paste("# Reading ", bedtools_output, sep = ""))
cov_per_bin <- read.delim(bedtools_output, stringsAsFactors = F, header = F)

# Two blacklists are created, one for the autosomes and one for chrX (because of the lower cov on chrX in males). These are merged in the final step
cov_per_bin_autosomal <- cov_per_bin[cov_per_bin[,1] %in% 1:29,]

# blacklist will contain the 2% bins with lowest coverage and 3% of bins with highest coverage
lcutoff <- quantile(cov_per_bin_autosomal[,4], 0.02)
ucutoff <- quantile(cov_per_bin_autosomal[,4], 0.97)

cov_per_bin_autosomal[,1] <- as.numeric(cov_per_bin_autosomal[,1])

# Plot the read count per bin with the cutoffs
ggplot(cov_per_bin_autosomal, aes(x = (V2+V3)/2, y = V4)) + geom_point(size = 0.5) + 
  facet_grid(~V1, space = "free", scales = "free") + 
  coord_cartesian(ylim = c(0, ucutoff*1.5)) + 
  geom_hline(aes(yintercept=lcutoff), color="red") +
  geom_hline(aes(yintercept=ucutoff), color="red") + theme_bw(base_size = 9) + 
  theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
  labs(x = "Genomic position", y = "Reads per bin")
ggsave(filename = paste(output_folder, blacklist_file, ".pdf", sep = ""), height = 3, width = 10)

cov_per_bin_X <- cov_per_bin[cov_per_bin[,1] =="X",]
lcutoff_X <- quantile(cov_per_bin_X[,4], 0.03)
ucutoff_X <- quantile(cov_per_bin_X[,4], 0.95)
ggplot(cov_per_bin_X, aes(x = (V2+V3)/2, y = V4)) + geom_point(size = 0.5) + 
  facet_grid(~V1, space = "free", scales = "free") + 
  coord_cartesian(ylim = c(0, ucutoff*2)) + 
  geom_hline(aes(yintercept=lcutoff_X), color="red") +
  geom_hline(aes(yintercept=ucutoff_X), color="red") + theme_bw(base_size = 9) + 
  theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
  labs(x = "Genomic position", y = "Reads per bin")
ggsave(filename = paste(output_folder, blacklist_file, "_X.pdf", sep = ""), height = 3, width = 10)

print(paste("# Writing blacklist ", output_folder, blacklist_file, sep = ""))
blacklist_auto <- cov_per_bin_autosomal[cov_per_bin_autosomal[,4] <= lcutoff | cov_per_bin_autosomal[,4] >= ucutoff,]
blacklist_X <- cov_per_bin_X[cov_per_bin_X[,4] <= lcutoff_X | cov_per_bin_X[,4] >= ucutoff_X,]
blacklist <- rbind(blacklist_auto, blacklist_X)

# Adjacent bins with a gap of one bin are merged together
blacklist_g <- GRanges(seqnames = blacklist[,1], IRanges(start = blacklist[,2]-as.numeric(bin_size)/2, end = blacklist[,3]+as.numeric(bin_size)/2))
blacklist_output <- reduce(blacklist_g)
start(blacklist_output) <- start(blacklist_output) + as.numeric(bin_size)/2
end(blacklist_output) <- end(blacklist_output) - as.numeric(bin_size)/2

sum(width(blacklist_output))

exportGRanges(blacklist_output, filename = paste(output_folder, blacklist_file, sep = ""), header = FALSE, chromosome.format = "NCBI")
