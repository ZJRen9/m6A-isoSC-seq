


def type_annotation_fileread(filename,fastadict,resultfile):
	"""
	chrom   position        snp_name        ref     alt     pvalue  ref_count       alt_count       baseSequence    alu_annotation exon_annotation  intron_annotation       editing_type    exon@intron     single@double   sequence_strand sequence_editing_type
	chr1    187418  .       G       A       35.77   94      15      GGGGGGG   False   False   False   intergentic     intergentic     intergentic     -      C-T
	chr1    944283  .       A       G       143.77  84      15      AAAAAAA     False   chr1;958458;944204;944204;944800;-;ENST00000477976#chr1;925730;944574;943907;944574;+;ENST00000342066   False   T-C     exon    double  -       T-C
	chr1    965490  .       C       T       124.77  23      8       cCCCCCC   chr1;960638;965602;964962;965602;+;ENST00000622660      False   C-T     exon    single  *       *
	"""
	f = open(filename,'r')
	d = open(resultfile,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			index_chorm = list_x.index("chrom")
			index_position = list_x.index("position")
			index_ref = list_x.index("ref")
			index_alt = list_x.index("alt")
			index_pvalue = list_x.index("pvalue")
			index_strand = list_x.index("sequence_strand")
			index_snptype = list_x.index("sequence_editing_type")
			str2write = str_x+"\t"+"motif"+"\n"
			d.write(str2write)
			continue
		chrom = list_x[index_chorm]
		position = list_x[index_position]
		ref = list_x[index_ref]
		alt = list_x[index_alt]
		pvalue = list_x[index_pvalue]
		strand = list_x[index_strand]
		snptype = list_x[index_snptype]
		if strand == "+":
			start = int(position) - 1
			end = start + 1
		if strand == "-":
			start = int(position) - 1
			end = start + 1
		if strand == "*":
			continue
		motif_start = start - 3
		motif_end = end + 3
		motif = fastadict[chrom][motif_start:motif_end]
		motif = sequence_reverse(sequence=motif,strand=strand)
		motif = motif[0:5]
		str2write = str_x+"\t"+motif + "\n"
		d.write(str2write)



def fasta_to_dict1(filename):
	"""
	>chr1
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	"""
	outdict = {}
	chromlist = []
	f = open(filename,'r')
	for str_x in f:
		str_x = str_x.strip("\n")
		if str_x[0] == ">":
			chrom = str_x[1:]
			chromlist.append(chrom)
			outdict[chrom] = []
		else:
			outdict[chrom].append(str_x)
			chromlist.append(chrom)
	print("step1 is over")
	for chrom in chromlist:
		sublist = outdict[chrom]
		substr = "".join(sublist)
		outdict[chrom] = substr
	print("step2 is over")
	return outdict

def fasta_to_dict2(filename):
	outdict = {}
	f = open(filename,'r')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		outdict[list_x[0]] = list_x[1]
	return outdict

def sequence_reverse(sequence,strand):
	reverse_dict = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
	base_list = []
	for i in range(len(sequence)):
		base = sequence[i]
		base = base.upper()
		if strand == "+":
			base_list.append(base)
		else:
			base = reverse_dict[base]
			base_list.append(base)
	if strand == "+":
		outstr = "".join(base_list)
	else:
		base_list.reverse()
		outstr = "".join(base_list)
	return outstr

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='motif add')
	parser.add_argument('--strand_annotation_result', required=True, help="strand_annotation_result")
	parser.add_argument('--fasta_file', required=True, help="fasta file")
	parser.add_argument('--result_file', required=True, help="merge result")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	strand_annotation_result = args.strand_annotation_result
	fasta_file = args.fasta_file
	result_file = args.result_file
	fasta_dict = fasta_to_dict2(filename=fasta_file)
	type_annotation_fileread(filename=strand_annotation_result,fastadict=fasta_dict,resultfile = result_file)

if __name__ == "__main__":
	main()

