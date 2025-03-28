#bin    swScore milliDiv    milliDel    milliIns    genoName    genoStart   genoEnd genoLeft    strand  repName repClass    repFa        mily   repStart    repEnd  repLeft id
#0   1892    83  59  14  chr1    67108753    67109046    -181847376  +   L1P5    LINE    L1  5301    5607    -544    1
#1   2582    27  0   23  chr1    8388315 8388618 -240567804  -   AluY    SINE    Alu -15 296 1   1
#1   4085    171 77  36  chr1    25165803    25166380    -223790042  +   L1MB5   LINE    L1  5567    6174    0   4
#1   2285    91  0   13  chr1    33554185    33554483    -215401939  -   AluSc   SINE    Alu -6  303 10  6
#1   2451    64  3   26  chr1    41942894    41943205    -207013217  -   AluY    SINE    Alu -7  304 1   8


def repeat_file_read(repeat_filename,genomesize,binsize):
	f = open(repeat_filename,'r')
	genomedict = empty_dictmake(genomesize,binsize)
	for str_x in f:
		if str_x[0]=="#":
			continue
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[5]
		chromstart = int(list_x[6])
		chromend = int(list_x[7])
		strand = list_x[9]
		repeat_name = list_x[10]
		repeat_class = list_x[11]
		repeat_family = list_x[12]
		bin_start = chromstart//binsize
		bin_end = chromend//binsize
		if repeat_family == "Alu":
			if bin_start == bin_end:
				try:
					genomedict[(chrom,bin_start)].append([chrom,chromstart,chromend,strand,repeat_name,repeat_class,repeat_family])
				except:
					pass
			else:
				for j in range(bin_start,bin_end+1):
					try:
						genomedict[(chrom,j)].append([chrom,chromstart,chromend,strand,repeat_name,repeat_class,repeat_family])
					except:
						pass
		else:
			pass
	return genomedict

def empty_dictmake(genomesize,binsize):
	'''
	chr1    248956422
	'''
	f = open(genomesize,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[0]
		totalsize = int(list_x[1])
		for i in range((totalsize//binsize)+1):
			outdict[chrom,i] = []
	return outdict

def counts_file_read(counts_filename,repeat_dict,binsize,annotation_result):
	'''
	chr1    16495   rs3210724       G       C       97.84   3       5
	chr1    51660   .       A       G       9.10    0       0
	chr1    51677   .       T       C       10.90   0       0
	'''
	f = open(counts_filename,'r')
	d = open(annotation_result,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0]=="chrom":
			d.write(str_x+"\t"+"alu_annotation"+"\n")
			chrom_index = list_x.index("chrom")
			position_index = list_x.index("position")
			continue
		chrom = list_x[chrom_index]
		pos = int(list_x[position_index])
		pos_bin = pos//binsize
		repeat_list = repeat_dict[chrom,pos_bin]
		overlap_repeat_list = repeat_overlap(pos,repeat_list)
		if overlap_repeat_list == []:
			Alu = "False"
		else:
			Alu_list = []
			for sub_overlap_repeat_list in overlap_repeat_list:
				subAlu = ";".join(sub_overlap_repeat_list)
				Alu_list.append(subAlu)
			Alu = "#".join(Alu_list)
		result_str = str_x+"\t"+Alu+"\n"
		d.write(result_str)

def repeat_overlap(pos,repeat_list):
	outlist = []
	for sublist in repeat_list:
		chrom,chromstart,chromend,strand,repeat_name,repeat_class,repeat_family = sublist
		if overlap(pos,chromstart,chromend)==True:
			outlist.append([chrom,str(chromstart),str(chromend),strand,repeat_name,repeat_class,repeat_family])
	return outlist

def overlap(pos,start,end):
	if pos >= start and pos <= end:
		return True
	else:
		return False

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--repeat_file', required=True, help="repeat file")
	parser.add_argument('--genome_size_file', required=True, help="genome size file")
	parser.add_argument('--fraction_file', required=True, help="fraction file")
	parser.add_argument('--annotation_result', required=True, help="annotation result")
	parser.add_argument('--binsize', required=False,default = "100000", help="binsize")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	repeat_filename = args.repeat_file
	genomesize = args.genome_size_file
	counts_filename = args.fraction_file
	annotation_result = args.annotation_result
	binsize = int(args.binsize)
	repeat_dict = repeat_file_read(repeat_filename,genomesize,binsize)
	counts_file_read(counts_filename,repeat_dict,binsize,annotation_result)


if __name__=="__main__":
	main()