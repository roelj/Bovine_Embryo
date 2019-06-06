# This script can be used to split BAM files into two strand-specific BAM files (FW/RV)
args <- commandArgs(trailingOnly=TRUE)
options(scipen = 999)

samplesheet_file <- args[1] # Requires columns: "Cell", "Run", "Method"
bam_folder <- args[2] # [bam_folder]/[run]/[cell].bam
output_folder <- args[3] # [output_folder]/FW/[run]/[cell].bam

# Read samplesheet
samplesheet <- read.delim(samplesheet_file, stringsAsFactors = F)

# Only select Strand-seq runs
strandseq_samples <- samplesheet[which(samplesheet$Method == "Strand-seq"),]

# Split the bam files using samtools view
for(i in 1:nrow(strandseq_samples)){
  
  bam_file <- paste(bam_folder, strandseq_samples[i,"Run"], "/", strandseq_samples[i,"Cell"], ".bam", sep = "")
  #print(bam_file)
  
  # Reads on FW strand can be selected using -F 20
  command_split_FW <- paste("samtools view -F 20 ",bam_file," -bh -o ",output_folder,"BAMs_FW_Strand/",strandseq_samples[i,"Run"], "/", strandseq_samples[i,"Cell"], ".bam", sep = "")
  
  # Reads on RV strand can be selected using -f 16
  command_split_RV <- paste("samtools view -f 16 ",bam_file," -bh -o ",output_folder,"BAMs_RV_Strand/",strandseq_samples[i,"Run"], "/", strandseq_samples[i,"Cell"], ".bam", sep = "")
  
  command_index_FW <- paste("samtools index ",output_folder,"BAMs_FW_Strand/",strandseq_samples[i,"Run"], "/",strandseq_samples[i,"Cell"], ".bam", sep = "")
  command_index_RV <- paste("samtools index ",output_folder,"BAMs_RV_Strand/",strandseq_samples[i,"Run"], "/",strandseq_samples[i,"Cell"], ".bam", sep = "")
  
  if(dir.exists(paste(output_folder,"BAMs_FW_Strand/",strandseq_samples[i,"Run"], sep = "")) == FALSE){
    dir.create(paste(output_folder,"BAMs_FW_Strand/",strandseq_samples[i,"Run"], sep = ""), recursive = T)
  }
  if(dir.exists(paste(output_folder,"BAMs_RV_Strand/",strandseq_samples[i,"Run"], sep = "")) == FALSE){
    dir.create(paste(output_folder,"BAMs_RV_Strand/",strandseq_samples[i,"Run"], sep = ""), recursive = T)
  }
  print(command_split_FW)
  system(command_split_FW)
  
  #print(command_index_FW)
  system(command_index_FW)
  
  print(command_split_RV)
  system(command_split_RV)
  
  #print(command_index_RV)
  system(command_index_RV)
}
print("# Done")
