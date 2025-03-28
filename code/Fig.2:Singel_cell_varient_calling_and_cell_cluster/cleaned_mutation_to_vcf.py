def bed_file_to_vcf(filename,outfilename):
	"""
	chr1    826948  826949  chr1_826948_826949_r    -       CCACT   22;5    NA;NA   NA;NA
	chr1    852099  852100  chr1_852099_852100_f    +       TTACA   77;7    NA;NA   NA;NA
	chr1    853904  853905  chr1_853904_853905_f    +       TGACC   48;6    54;7    NA;NA
	chr1    944080  944081  chr1_944080_944081_f    +       AAACT   50;10   NA;NA   NA;NA
	chr1    944207  944208  chr1_944207_944208_f    +       CAACA   43;9    NA;NA   NA;NA
	chr1    944209  944210  chr1_944209_944210_f    +       ACACA   46;7    NA;NA   NA;NA
	chr1    944579  944580  chr1_944579_944580_r    -       ACACG   91;11   NA;NA   103;13
	"""
	f = open(filename,'r')
	d = open(outfilename,'a')
	outdict = {}
	outlist = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[0]
		start = list_x[1]
		end = list_x[2]
		name = list_x[3]
		strand = list_x[4]
		if strand == "+":
			ref = "C"
			alt = "T"
		else:
			ref = "G"
			alt = "A"
		#############################################################
		str2write = "\t".join([chrom,start,".",ref,alt,"100",".","AC=1;AF=0.500;AN=2;BaseQRankSum=-5.668;DP=71;Dels=0.00;ExcessHet=3.0103;FS=1.203;HaplotypeScore=36.1509;MLEAC=1;MLEAF=0.500;MQ=60.00;MQ0=0;MQRankSum=0.000;QD=0.95;ReadPosRankSum=1.113;SOR=0.489","GT:AD:DP:GQ:PL","0/1:59,12:71:96:96,0,1633"]) + "\n"
		d.write(str2write)
		#############################################################
	d.close()

def header_write(filename,outfilename):
	f = open(filename,'r')
	d = open(outfilename,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if str_x[0] == "#":
			str2write = str_x + "\n"
			outfile.write(str2write)
		else:
			pass
	d.close()

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='convert mutation to vcf')
	parser.add_argument("--vcf_header_filename", required = False, help = "vcf file path")
	parser.add_argument("--bed_filename", required = False, help = "mutation file path")
	parser.add_argument("--output_vcf_filename", required = False, help = "mutation file path")
	args = parser.parse_args()
	return args


if __name__ == "__main__":
	args = make_args()
	bed_filename = args.bed_filename
	vcf_filename = args.vcf_filename
	output_vcf_filename = args.output_vcf_filename
	header_write(filename=vcf_filename,outfilename=output_vcf_filename)
	bed_file_to_vcf(filename=bed_filename,outfilename=output_vcf_filename)
	#resultwrite(outlist,outdict)
