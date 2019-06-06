# This script selects all heterozygous SNPs in the input VCF.
# filters all SNPs from a VCF that have less than [-c N] reads on the alternative or reference allele
# It also removes indels and positions that have a depth of more than [-d N].

import getopt
import sys
import re
import gzip
import os

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:c:z:d:', ['input=','output=','zip=', 'counts=','depth='])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)

	elif opt in ('-i', '--input'):
		input_vcf = arg
	elif opt in ('-o', '--output'):
		output_file = arg
	elif opt in ('-c', '--counts'):
		min_counts = arg
	elif opt in ('-z', '--zip'):
		zip_output = arg
	elif opt in ('-d', '--zip'):
		max_depth = arg

if re.search(".gz", input_vcf):		
	f1 = gzip.open(input_vcf)
else:
	f1 = open(input_vcf)

# Only variants on chr1 to chr29 and chrX are included
chromosomes = map('{:1}'.format, range(1,30))
chromosomes.append("X")

print "# Selecting heterozygous SNPs from " + input_vcf
print "# Writing output to " + output_file

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
				DP = map(int, re.match(".*DP=-*(\d+);", columns[7]).group(1).split(","))
				DP4 = map(int, re.match(".*DP4=-*(\d+,\d+,\d+,\d+);", columns[7]).group(1).split(","))
				REF_counts = DP4[0]+DP4[1]
				ALT_counts = DP4[2]+DP4[3]
				if int(REF_counts) > int(min_counts) and int(ALT_counts) > int(min_counts) and int(DP[0]) < int(max_depth):
					#print line + "\n"
					#print REF_counts
					#print ALT_counts
					f2.write(line + "\n")
	else:
		# write headers to output file
		line = line.rstrip()
		f2.write(line + "\n")

print "# Filtering Done"

f1.close()
f2.close()

if zip_output == "TRUE":
	print "# Zipping output file"
	print "bgzip -c " + output_file + " > " + output_file + ".gz"
	os.system("bgzip -c " + output_file + " > " + output_file + ".gz")
	print "bcftools index " + output_file + ".gz"
	os.system("bcftools index " + output_file + ".gz")
	print "bcftools stats " + output_file + ".gz"
	os.system("bcftools stats " + output_file + ".gz > " + output_file + ".vchk")

