# This script is used to annotate multiple VCFs from single cells. 
# The annotation file should be sorted, bgzipped and indexed (tabix). It should contain the columns: chr-start-end-annotation
# The info_field will be the name of the annotation

args <- commandArgs(trailingOnly=TRUE)
options(scipen = 999)

samplesheet_file <- args[1]
vcf_folder <- args[2] # directory containing input VCFs [vcf_folder]/[Run]/*.vcf
annotation_file <- args[3] # File should contain chr-start-end-info-to-add
header_file <- args[4] # txt file containing line that will be added to the VCF header, eg: ##INFO=<ID=PAT,Number=1,Type=Integer,Description="Genotype of father">
output_folder <- args[5]
info_field <- args[6] # this ID will be added to the INFO field in the output VCF
runs <- args[7] # Leave empty to run for all sequencing runs

dir.create(output_folder)

# Unless a specific run is specified, VCF for all samples listed in the samplesheet will be annotated (if the VCF exists)
samplesheet <- read.delim(samplesheet_file, stringsAsFactors = F)
print(runs)

if(!is.na(runs)){
  samplesheet <- samplesheet[which(samplesheet$Run %in% runs),]
}

## Annotate single cell VCFs
if(file.exists(annotation_file)){
  if(file.exists(header_file)){
    for(Cell in unique(samplesheet$Cell)){
      
      VCF_file <- list.files(vcf_folder, pattern = "\\.vcf.gz$", recursive = T)[grep(pattern = Cell, x = list.files(vcf_folder, pattern = "\\.vcf.gz$", recursive = T))]
      
      # Only run if VCF file exists:
      if(file.exists(paste(vcf_folder, VCF_file, sep = "")) == TRUE & file_test("-f", paste(vcf_folder, VCF_file, sep = ""))){
        output_folder_run <- paste(output_folder,samplesheet$Run[samplesheet$Cell == Cell], "/", sep = "")
        
        if(dir.exists(output_folder_run) == FALSE){
          dir.create(output_folder_run, showWarnings = F, recursive = T)
        }
        
        annotation_command <- paste("bcftools annotate -a ",annotation_file," -c CHROM,FROM,TO,INFO/",info_field," -h ",header_file," -O z -o ", output_folder_run, Cell,".vcf.gz  ",vcf_folder, VCF_file, sep = "")
        print(annotation_command)
        system(annotation_command)
      } else {
          print(paste("! VCF File not found: ", VCF_file, sep = ""))
        }
    }
  } else {
    print(paste("! Header file not found: ", header_file, sep = ""))
    
  }

} else {
  print(paste("! Annotation file not found: ", annotation_file, sep = ""))
}

