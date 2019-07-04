# This script will run idxstats to determine the number of reads on each chromosome for each sample
args <- commandArgs(trailingOnly=TRUE)
options(scipen = 999)

# The samplesheet should contain the following columns: Run (eg the name for the folder containing the bams), Method (Strand-seq or WGS)
samplesheet_file <- args[1]
# Input folder should contain subfolders for each run containing bam files. [bam_folder]/[run]/*.bam
input_folder <- args[2] 
# idxstats data will be written here [output_folder]/idxstats/run/*.idxstats
output_folder <- args[3]
# The script will only run for one folder if "runs" is defined (use name folder containing bam files)
runs <- args[4]

samplesheet <- read.delim(samplesheet_file, header = T, stringsAsFactors = F)
# Remove positive and negative controls:
samplesheet <- samplesheet[samplesheet$Embryo != "POS" & samplesheet$Embryo != "NEG",]

# only select the specified runs
if(!is.na(runs)){
  samplesheet <- samplesheet[samplesheet$Run %in% runs,]
}

dir.create(output_folder, showWarnings = F)

print("## Running Samtools idxstats")
for(run in unique(samplesheet$Run)){
  for(bam_file in list.files(paste(input_folder, run, sep = ""), pattern = "\\.bam$")){
    input_file <- paste(input_folder, run, "/", bam_file, sep ="")
    cell <- gsub(x = bam_file, pattern = ".bam", replacement = "")
    idxstats_folder <- paste(output_folder, "idxstats/", run, "/", sep = "")
    if(dir.exists(idxstats_folder) == FALSE){
      dir.create(idxstats_folder, showWarnings = F, recursive = T)
    }
    
    if(file.exists(paste(idxstats_folder, cell, ".idxstats", sep  ="")) == FALSE | file.info(paste(idxstats_folder, cell, ".idxstats", sep  =""))$size == 0){
      command <- paste("samtools idxstats ", input_file, " > ", idxstats_folder, cell, ".idxstats", sep  ="")
      print(command)
      system(command)
    }
  }
}