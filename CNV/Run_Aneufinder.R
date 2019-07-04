# This script is used to run Aneufinder on all bam files in one folder.
# Requires AneuFinder, BSgenome.Btaurus.UCSC.bosTau8 and/or BSgenome.Hsapiens.UCSC.hg19 packages
# Can be run by run_aneufinder.sh script

args <- commandArgs(trailingOnly=TRUE)
bamfile_folder <- args[1]
output_folder <- args[2]
bin_size_input <- args[3]
strand_seq <- args[4]
reference <- args[5]
blacklist <- args[6]
methods_input <- args[7]
cnv_states_input <- args[8]
chromosome_input <- args[9]

library(AneuFinder)

if(reference == "BSgenome.Btaurus.UCSC.bosTau8"){
  library(BSgenome.Btaurus.UCSC.bosTau8)
  genome <- BSgenome.Btaurus.UCSC.bosTau8
  if(is.na(chromosome_input)){
    chromosomes <- c(1:29, "X")
  }
} else if(reference == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  if(is.na(chromosome_input)){
    chromosomes <- c(1:22, "X")
  }
}

# need to check if chromosomes need to be numeric or can be string (else "X" can be problematic?)
if(!is.na(chromosome_input)){
  chromosomes <- unlist(strsplit(chromosome_input, split = ","))
}

print(paste("### Running aneufinder for ", length(list.files(bamfile_folder, pattern = "\\.bam$")), " BAM files", sep = ""))

print(paste("### Input folder = ", bamfile_folder, sep = ""))

print(paste("### Output folder = ", output_folder, sep = ""))

print(paste("### Strand-seq = ", strand_seq, sep = ""))


if(!is.na(cnv_states_input)){
  cnv.states <- unlist(strsplit(cnv_states_input, split = ","))
} else {
  cnv.states <- c("zero-inflation", paste0(0:10, "-somy"))
}

if(!is.na(methods_input)){
  methods <- unlist(strsplit(methods_input, split = ","))
} else {
  methods <- c("HMM", "dnacopy", "edivisive")
}


print(paste("### cnv.states = ", paste(cnv.states, collapse = ","), sep = ""))

bin_size <- unlist(strsplit(bin_size_input, split = ","))

dir.create(output_folder)

Aneufinder(inputfolder = bamfile_folder, outputfolder = output_folder, assembly = NULL,
           numCPU = 4, states = cnv.states, chromosomes = chromosomes, remove.duplicate.reads = T, 
           correction.method = c("GC"), binsizes = as.numeric(bin_size), reads.store=FALSE, 
           blacklist = blacklist, strandseq = as.logical(strand_seq),
           GC.BSgenome = genome, method = methods, stepsize = 4e4, min.mapq = 10)

print("### Finished Aneufinder analysis")

