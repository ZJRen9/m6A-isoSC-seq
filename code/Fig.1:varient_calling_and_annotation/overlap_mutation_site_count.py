import argparse

def overlap_fileread(filename):
	"""
	chr1    29553   29703   ENSG00000243485_MIR1302-2HG_0_0 +       .       -1      -1      .       .       .       .       .  0
	chr1    29563   29713   ENSG00000243485_MIR1302-2HG_0_1 +       .       -1      -1      .       .       .       .       .  0
	chr1    29573   29723   ENSG00000243485_MIR1302-2HG_0_2 +       .       -1      -1      .       .       .       .       .  0
	chr1    1254738 1254888 ENSG00000160087_UBE2J2_0_83 -   chr1    1254887 1254888 chr1_1254887_1254888_r  TCACT   -            37  24  1
	chr1    1254748 1254898 ENSG00000160087_UBE2J2_0_84 -   chr1    1254887 1254888 chr1_1254887_1254888_r  TCACT   -            37  24  1
	"""
	f = open(filename,'r')
	outdict1 = {}
	outdict2 = {}
	outdict3 = {}
	outlist = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[0]
		start = list_x[1]
		end = list_x[2]
		windowname = list_x[3]
		strand = list_x[5]
		overlen = int(list_x[-1])
		###################################
		try:
			outdict1[windowname] += 0
			outdict2[windowname] += 0
			outdict3[windowname] += 0
		except:
			outdict1[windowname] = 0
			outdict2[windowname] = 0
			outdict3[windowname] = 0
			outlist.append([chrom,start,end,windowname,strand])
		###################################
		if overlen >= 1:
			outdict1[windowname] +=1
			ref_count = int(list_x[-3])
			alt_count = int(list_x[-2])
			outdict2[windowname] += ref_count
			outdict3[windowname] += alt_count
		else:
			outdict1[windowname] +=0
			outdict2[windowname] +=0
			outdict3[windowname] +=0
		###################################
	return outdict1,outdict2,outdict3,outlist

def resultwrite(outdict1,outdict2,outdict3,outlist,output_windown_mutation_count):
	d = open(output_windown_mutation_count,'a')
	for sublist in outlist:
		chrom,start,end,windowname,strand = sublist
		value1 = outdict1[windowname]
		value2 = outdict2[windowname]
		value3 = outdict3[windowname]
		str2write = "\t".join([chrom,start,end,windowname,strand,str(value1),str(value2),str(value3)]) + "\n"
		d.write(str2write)

def make_args():
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-I', '--input_intersect_result', required=True, help="input intersect result")
	parser.add_argument('-O', '--output_windown_mutation_count', required=True, help="output windown mutation count")
	args = parser.parse_args()
	return args


def main():
	args = make_args()
	input_intersect_result = args.input_intersect_result
	output_windown_mutation_count = args.output_windown_mutation_count 
	########################################
	outdict1,outdict2,outdict3,outlist = overlap_fileread(filename=input_intersect_result)
	resultwrite(outdict1=outdict1,outdict2=outdict2,outdict3=outdict3,outlist=outlist,output_windown_mutation_count=output_windown_mutation_count)

if __name__=="__main__":
	main()


