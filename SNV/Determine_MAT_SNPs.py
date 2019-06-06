# This script collects all PAT_MAT and MAT SNPs from a VCF file and writes the positions to a bed file

import getopt
import sys
import re
import gzip
import os
import fnmatch

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:z:', ['input=','output=','zip='])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)

	elif opt in ('-i', '--input'):
		input_folder = arg
	elif opt in ('-o', '--output'):
		output_folder = arg

# Add a slash if the output folder path does not contain one yet
if output_folder.endswith("/") == False:
	output_folder = output_folder + "/"

# Loop to all VCF files in the input folder (recursively)
VCF_files = []
for root, dirnames, filenames in os.walk(input_folder):
    for filename in fnmatch.filter(filenames, '*.vcf.gz'):
        VCF_files.append(os.path.join(root, filename))

for input_vcf in VCF_files:
	print input_vcf
	if re.search(".gz", input_vcf):		
		f1 = gzip.open(input_vcf)
	else:
		f1 = open(input_vcf)

	# Only variants on chr1 to chr29 and chrX are included
	chromosomes = map('{:1}'.format, range(1,30))
	chromosomes.append("X")
	
	output_name = re.sub(input_folder, "", input_vcf)
	output_name = re.sub(".vcf.gz", ".bed", output_name)
	
	run_temp = input_vcf.split("/")
	run = run_temp[len(run_temp)-2]
	output_folder_file = output_folder + run + "/"
	if not os.path.exists(output_folder_file):
		os.makedirs(output_folder_file)

	output_file = output_folder + output_name
	
	print output_file

	# Output will be written to f2
	f2 = open(output_file, "w")

	for line in f1:
		if not line.startswith("#"):
			line = line.rstrip()
			columns = line.split("\t")
			chromosome = columns[0]
			if chromosome in chromosomes:
				info = columns[7]
				# Remove indels:
				if not "INDEL" in info:
					if "PAT=" in info:
						PAT = map(int, re.match(".*PAT=-*(\d+,\d+)", columns[7]).group(1).split(","))
						GT = columns[9].split(":")
									
						output = [chromosome, columns[1],columns[1]]
						#print output
						#print line
						#print PAT
						#print GT[0]
						if PAT[0] == 0 and PAT[1] == 0 and GT[0] == "0/0":
							#print	"PAT/MAT"
							output.append("PAT/MAT")
						elif PAT[0] == 1 and PAT[1] == 1 and GT[0] == "1/1":
							#print	"PAT/MAT"
							output.append("PAT/MAT")
						else:
							#print "MAT"
							output.append("MAT")
						
						output2 = "\t".join(output)
						#print output2
						f2.write(output2 + "\n")

			
		else:
			# write headers to output file
			line = line.rstrip()
			#f2.write(line + "\n")

	f1.close()
	f2.close()

print "# Filtering Done"
