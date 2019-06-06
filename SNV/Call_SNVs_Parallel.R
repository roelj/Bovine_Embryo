# This script will generate jobs to call SNVs for multiple sequencing runs
args <- commandArgs(trailingOnly=TRUE)
options(scipen = 999)

# SNVs can be genotyped for SNP positions described by a regions file, or SNVs different from reference can be called genome wide 
# Options: TRUE / FALSE
Genotype <- args[1]

# The samplesheet should contain the following columns: Run (eg the name for the folder containing the bams), Method (Strand-seq or WGS)
samplesheet_file <- args[2]
# Input folder should contain the Call_SNVs.R script
input_folder <- args[3] 
# The bam_folder should contain subfolders for each run containing bam files. [bam_folder]/[run]/*.bam
bam_folder <- args[4]
# SNV VCFs will be written to this folder:
output_folder <- args[5] # [output_folder][run]/[library].vcf.gz
# regions_file should be a bgzipped, sorted tab-sep .txt file containing columns chr-start-end
regions_file <- args[6]

# path to the refence genome (Bos_taurus.UMD3.1.dna.toplevel.fa)
reference_path <- args[7]

# Jobs/logs will be called after job_name (eg Call_SNPs[_Run_Part])
job_name <- args[8]

# Will only be run for the run specified by 7th argument. Leave empty to run for all samples
runs <- args[9]

if(is.na(job_name)){
  job_name <- "Call_SNVs"
}

## Settings
# Jobs and logs will be written to these folders:
logs_folder <- paste(input_folder, "Logs/", sep = "")
jobs_folder <- paste(input_folder, "Jobs/", sep = "")
dir.create(logs_folder, showWarnings = F, recursive = T)
dir.create(jobs_folder, showWarnings = F, recursive = T)

samplesheet <- read.delim(samplesheet_file, header = T, stringsAsFactors = F)
# Remove positive and negative controls:
samplesheet <- samplesheet[samplesheet$Embryo != "POS" & samplesheet$Embryo != "NEG",]

# only select the specified runs
if(!is.na(runs)){
  samplesheet <- samplesheet[samplesheet$Run %in% runs,]
}

dir.create(output_folder, showWarnings = F, recursive = T)

Call_SNV_Script <- paste(input_folder, "Call_SNVs.R", sep = "")
print(Call_SNV_Script)
if(file.exists(Call_SNV_Script)){
  if(file.exists(reference_path)){
    print(unique(samplesheet$Run))
    for(run in unique(samplesheet$Run)){
      run_info <- samplesheet[which(samplesheet$Run == run),]
      bam_folder_run <- paste(bam_folder, run, "/", sep = "")
      
      if(dir.exists(bam_folder_run) == TRUE){
        
      samples_run <- samplesheet[which(samplesheet$Run == run),]
      batches <- ceiling(nrow(samples_run)/12)
      
      # SNPs will be called per batch of 12 libraries to speed up the process (one job per 12 libraries)
      for(i in 0:(batches-1)){
        print(i)
        
        if(i*12+12 < nrow(samples_run)){
          samples <- paste(samples_run$Cell[(i*12+1) : (i*12+12)], collapse = ",")
          print(samples)
        } else {
          samples <- paste(samples_run$Cell[(i*12+1) : nrow(samples_run)], collapse = ",")
          print(samples)
        }
      
       print(bam_folder_run)
       if(length(list.files(bam_folder_run, pattern = "\\.bam$")) > 0){
         job_ID <- paste(jobs_folder, job_name, "_", run, "_part", (i+1) , ".sh", sep = "")
         shell_script <- paste("#!/bin/bash

#$ -S /bin/bash
#$ -l h_vmem=30G
#$ -l h_rt=03:00:00
#$ -cwd
#$ -o ",logs_folder, format(Sys.Date(), "%Y%m%d"),"_",job_name,"_", run, "_part", (i+1),"_log.txt
#$ -e ",logs_folder, format(Sys.Date(), "%Y%m%d"),"_",job_name,"_", run, "_part", (i+1), "_log.txt

GENOTYPE=",Genotype,"

# Folder containing all the bam files (subfolder should start with the Run_ID eg BAM_FOLDER/RUN/{}.bam
BAM_FOLDER=",bam_folder,"

# Output will be written to this folder:
OUTPUT_FOLDER=",output_folder,"

# Path to the samplesheet. Should contain ID (= cell id), Run
SAMPLESHEET=",samplesheet_file,"

# Script is ran by sequencing run:
RUN=",run,"

# 
MIN_QUAL=10

# Path to reference
REFERENCE_PATH=",reference_path,"

REGIONS=",regions_file,"

SAMPLES=\"",samples,"\"

guixr load-profile ~/.guix-profile/ -- <<EOF

	Rscript ",Call_SNV_Script," $GENOTYPE $BAM_FOLDER $OUTPUT_FOLDER $SAMPLESHEET $RUN $MIN_QUAL $REFERENCE_PATH $REGIONS $SAMPLES

EOF
", sep = "")
        
        print(paste("# Writing shell script: ",job_ID, sep = ""))
        write.table(shell_script, file = job_ID, sep = "\t", quote = F, row.names = F, col.names = F)
        
        command <- paste("qsub ", job_ID, sep ="")
        print(command)
        system(command)
    } else {
      print("! No BAM files found in folder: ",bam_folder_run, " !", sep = "")
    }
    }
      } else {
        print(paste("! Folder not found: ", bam_folder_run, " !"), sep = "")
        print(paste("! Skipping run ", run, sep = ""))
      }
    }
    
    } else {
      print(paste("! Reference genome not found: ", reference_path, " !"), sep = "")
    }
  } else {
    print(paste("! SNV calling script: ", Call_SNV_Script, " not found!"), sep = "")
  }


  

