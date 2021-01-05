# This script can be used to parallize genotyping by splitting the regions file by chromosome

args <- commandArgs(trailingOnly=TRUE)
input_VCF <- args[1] # full path to unzipped vcf
input_BAM <- args[2] # full path to bam file
input_folder <- args[3] # full path to folder containing scripts, jobs and logs will be written here: [input_folder]/Jobs/ and [input_folder]/Logs
output_folder <- args[4] # full path to output folder
output_name <- args[5] # name of output file [output_folder][output_name].vcf
reference <-  args[6] # full path to reference genome

dir.create(output_folder)

if(sub('.*(?=.$)', '', x=input_folder, perl=T) != "/"){
  input_folder <- paste(input_folder, "/", sep = "")
}

jobs_folder <- paste(input_folder, "Jobs/", output_name, "/", sep = "")
dir.create(jobs_folder, showWarnings = F, recursive = T)
logs_folder <- paste(input_folder, "Logs/", output_name, "/", sep = "")
dir.create(logs_folder, showWarnings = F, recursive = T)

# The temp_folder will not be automatically removed
temp_folder <- paste(output_folder, "temp/", sep = "")

dir.create(temp_folder, recursive = T)

## Obtain the coordinates of the SNP positions (chr:start) from the VCF
print("# Obtain SNP positions from input VCF")
vcf_to_positions <- paste("sed -e 's/chr//' ",input_VCF," | awk '{OFS=\"\t\"; if (!/^#/){print $1,$2}}' > ",output_folder,"SNP_Positions.txt", sep = "")
print(vcf_to_positions)
system(vcf_to_positions)

## Split the SNP positions to a seperate file per chr
print("# Split SNP positions per chromosome")
split_vcf_command <- paste("awk '{print $0 >> \"",temp_folder,"\"$1\".txt\"}' ",output_folder,"SNP_Positions.txt", sep = "")
print(split_vcf_command)
system(split_vcf_command)

print("# Generating jobs")

for(chr_file in list.files(temp_folder)){
  
  chr <- gsub(".txt", replacement = "", x= chr_file)
  job_ID <- paste(jobs_folder, "Call_SNPs_chr", chr, ".sh", sep = "")
  
  output_file <- paste(temp_folder, "SNPs_chr", chr, ".vcf.gz", sep = "")
  shell_script <- paste("
#!/bin/bash

COMMAND=\"bcftools mpileup -Ou -f ",reference," -q 10 -Q 20 -d 10000 -R ",temp_folder, chr_file," ",input_BAM," | bcftools call -mO z -P 0 -o ",output_file,"\"
echo \"$COMMAND\"
eval \"$COMMAND\"
echo \"Finished\"
", sep = "")
  
  print(paste("# Writing shell script: ",job_ID, sep = ""))
  write.table(shell_script, file = job_ID, sep = "\t", quote = F, row.names = F, col.names = F)
  
  command <- paste("qsub -V -l h_vmem=30G -l h_rt=04:00:00 -cwd -o ",logs_folder, format(Sys.Date(), "%Y%m%d"),"_Call_SNPs_chr", chr,"_log.txt -e ",logs_folder, format(Sys.Date(), "%Y%m%d"),"_Call_SNPs_chr", chr,"_log.txt ", job_ID, sep ="")
  print(command)
  system(command)
}

print("# Running jobs")


