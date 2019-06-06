## This script is used to call single nucleotide variants in the single cell sequencing data.
# Can be ran via Call_SNVs.sh

args <- commandArgs(trailingOnly=TRUE)
genotype <- args[1] # TRUE/FALSE
bam_folder <- args[2] # path to directory containing subfolders per run containing the BAM files [bam_folder]/[run]/*.bam
output_folder <- args[3]
samplesheet_file <- args[4] #Samplesheet should contain Library IDs and Run IDs
run <- args[5]
min_qual <- args[6] # Minimal mapping quality of reads included for variant calling
reference <- args[7] # full path to reference genome (.fa)
regions_file <- args[8] # full path to tab seperated file (pref sorted, bgzipped and tabixed)
samplelist <- args[9] # Cell ids seperated by comma

min_base_quality <- 20 # minimal base quality used to call SNVs (bcftools default is 13, but generates false positives)

overwrite = FALSE

print(paste("## Calling SNVs for run: ",run, sep = ""))
dir.create(output_folder)
dir.create(paste(output_folder, run, "/", sep = ""), showWarnings = F)
stats_folder <- paste(output_folder, "VCF_Stats/", sep = "")
dir.create(stats_folder)

samplesheet <- read.delim(samplesheet_file, stringsAsFactors = F)
samples_run <- samplesheet$Library[samplesheet$Run == run]
print(samplelist)
if(samplelist != ""){
  sample_IDs <- unlist(strsplit(samplelist, split = ","))
  samples <- samplesheet$Cell[which(samplesheet$Cell %in% sample_IDs)]
  print(samples)
} else {
  samples <- samples_run
}

## SNP calling whole genome
if(genotype == "FALSE"){
  print("## Calling SNVS")
  for(cell in samples){
    #run <- samplesheet$Run[samplesheet$Library == cell]
    
    output_file <- paste(output_folder, run, "/", cell, ".vcf.gz", sep = "")
    
    dir.create(paste(output_folder, run, "/", sep = ""), showWarnings = F)
    
    command <- paste("bcftools mpileup -Ou -f ",reference," -q ", min_qual, " ",bam_folder, run, "/", cell, ".bam | bcftools call -vmO z -P 0 -o ",output_file, sep = "")
    print(command)
    system(command)
    
    idx_command <- paste("bcftools index ", output_file, sep = "")
    print(idx_command)
    system(idx_command)
    
    stats_command <- paste("bcftools stats ",output_file," > ",stats_folder, cell,".vchk" ,sep = "")
    print(stats_command)
    system(stats_command)
  }
}

## Genotype ech single cell across positions in the regions_file
if(genotype == "TRUE"){
  
  regions <- regions_file 
  print(paste("## Genotyping for: ", regions, sep = ""))

  if(file.exists(regions)){
    for(cell in samples){
      print(paste("# ", cell, sep = ""))

      output_file <- paste(output_folder, run, "/", cell, ".vcf.gz", sep = "")
      
      if(overwrite == FALSE | file.exists(output_file) == FALSE){
        command <- paste("bcftools mpileup -Ou -f ",reference," -q ", min_qual, " -Q ",min_base_quality," -R ",regions, " ", bam_folder, run, "/", cell, ".bam | bcftools call -mO z -P 0 -V indels -o ",output_file, sep = "")
        print(command)
        system(command)
      
        idx_command <- paste("bcftools index ", output_file, sep = "")
        print(idx_command)
        system(idx_command)
        
        if(overwrite == FALSE | file.exists(paste(stats_folder, cell,".vchk", sep = "")) == FALSE){
          stats_folder_run <- paste(stats_folder, run, "/", sep = "")
          if(dir.exists(stats_folder_run) == FALSE){
            dir.create(stats_folder_run)
          }
          stats_command <- paste("bcftools stats ",output_file," > ",stats_folder_run, cell,".vchk" ,sep = "")
          print(stats_command)
          system(stats_command)
        }
      }
      }
    } else {
      print(paste("! Regions file not found: ", regions, sep = ""))
    }
}

print("### Finished SNV calling")

