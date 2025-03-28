

def gene_annotation_file_read(annotation_file,result_file,sequence_type):
	"""
	chr1    8878693 .       G       C       326.77  513     3       False   chr1;8878720;8868026;8878579;8878720;-;ENST00000646156  chr1;8879190;8868048;8874917;8878811;-;ENST00000489867   C-G     exon    single
	chr1    8878696 .       G       C       67.77   459     1       False   chr1;8878720;8868026;8878579;8878720;-;ENST00000646156  chr1;8879190;8868048;8874917;8878811;-;ENST00000489867   C-G     exon    single
	chr1    8878709 .       T       C       12.99   335     0       False   chr1;8878720;8868026;8878579;8878720;-;ENST00000646156  chr1;8879190;8868048;8874917;8878811;-;ENST00000489867   A-G     exon    single
	"""
	f = open(annotation_file,'r')
	d = open(result_file,'a')
	#header_str = "\t".join(["chrom","position","site","ref","alt","pvalue","ref_depth","alt_depth","alu_annotation","exon_annotation","intron_annotation","editing_type","exon@intron","single@double"]) + "\n"
	#d.write(header_str)
	reverse_dict = {"A":"T","T":"A","C":"G","G":"C"}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			#print(list_x)
			chrom_index = list_x.index("chrom")
			position_index = list_x.index("position")
			#site_index = list_x.index("site")
			site_index = list_x.index("snp_name")
			ref_index = list_x.index("ref")
			alt_index = list_x.index("alt")
			baseSequence_index = list_x.index("baseSequence")
			alu_index = list_x.index("alu_annotation")
			exon_index = list_x.index("exon_annotation")
			intron_index = list_x.index("intron_annotation")
			d.write("\t".join([str_x,"editing_type","exon@intron","single@double","sequence_strand","sequence_editing_type"])+"\n")
			continue
		chrom = list_x[chrom_index]
		site_name = list_x[position_index]
		site = list_x[site_index]
		ref = list_x[ref_index]
		alt = list_x[alt_index]
		if len(ref) > 1 or len(alt) > 1:
			continue
		baseSequence = list_x[baseSequence_index]
		overlap_alu_str = list_x[alu_index]
		overlap_exon_str = list_x[exon_index]
		overlap_intron_str = list_x[intron_index]

		exon_strand,exon_strand_type = exon_strand_estimate(overlap_exon_str)
		intron_strand,intron_strand_type = intron_strand_estimate(overlap_intron_str)
		alu_strand,alu_strand_type = alu_strand_estimate(overlap_alu_str)
		sequence_strand = baseSequence_strand_estimate(baseSequence,sequence_type)
		#print(overlap_exon_str)
		#print(exon_strand,exon_strand_type)
		#print(overlap_intron_str)
		#print(intron_strand,intron_strand_type)
		strand_type,snp_strand,local = snp_strand_estimate(exon_strand,exon_strand_type,intron_strand,intron_strand_type)
		if snp_strand == "+":
			snp_type = ref + "-" + alt
			result_str = str_x + "\t" +"\t".join([snp_type,local,strand_type])
		elif snp_strand == "-":
			snp_type = reverse_dict[ref] + "-" + reverse_dict[alt]
			result_str = str_x + "\t" +"\t".join([snp_type,local,strand_type])
		elif snp_strand == "intergentic":
			#snp_type = ref + "-" + alt + ";" + reverse_dict[ref] + "-" + reverse_dict[alt]
			snp_type = "intergentic"
			result_str = str_x + "\t" +"\t".join([snp_type,local,strand_type])
		else:
			raise Exception("ERROE")

		if sequence_strand == "+":
			snp_type = ref + "-" + alt
			result_str = result_str + "\t" +"\t".join([sequence_strand,snp_type]) + "\n"
		elif sequence_strand == "-":
			snp_type = reverse_dict[ref] + "-" + reverse_dict[alt]
			result_str = result_str + "\t" +"\t".join([sequence_strand,snp_type]) + "\n"
		elif sequence_strand == "*":
			snp_type = "*"
			result_str = result_str + "\t" +"\t".join([sequence_strand,snp_type]) + "\n"
		else:
			raise Exception("ERROE")
		d.write(result_str)

def baseSequence_strand_estimate(baseSequence,sequence_type):
	"""
	aaaaggaaaaagagggaaagaaagaagggaaaagaaaaaggaggagaaaa
	TTTTTTTTTTTTTTtCTTTTTTTTTTTTTTTTTTTTTTTCTTCCTCTTCCTCTTTTTTT
	"""
	forword_counts = 0
	reverse_counts = 0
	for i in range(len(baseSequence)):
		base_i = baseSequence[i]
		if base_i in ['a','t','c','g']:
			reverse_counts += 1
		elif base_i in ['A','T','C','G']:
			forword_counts += 1
		else:
			pass
	total_counts = forword_counts + reverse_counts
	#print(total_counts,forword_counts,reverse_counts)
	if sequence_type == "fr-strand":
		if forword_counts > 0.7*total_counts and forword_counts <= 1.0*total_counts:
			strand_type = "+"
		elif forword_counts > 0.3*total_counts and forword_counts <= 0.7*total_counts:
			strand_type = "*"
		elif forword_counts >= 0.0*total_counts and forword_counts <= 0.3*total_counts:
			strand_type = "-"
		else:
			raise Exception("ERROE")
	elif sequence_type == "rf-strand":
		if forword_counts > 0.7*total_counts and forword_counts <= 1.0*total_counts:
			strand_type = "-"
		elif forword_counts > 0.3*total_counts and forword_counts <= 0.7*total_counts:
			strand_type = "*"
		elif forword_counts >= 0.0*total_counts and forword_counts <= 0.3*total_counts:
			strand_type = "+"
		else:
			raise Exception("ERROE")
	else:
		strand_type = "*"
	return strand_type




def snp_strand_estimate(exon_strand,exon_strand_type,intron_strand,intron_strand_type):
	#print([exon_strand,exon_strand_type,intron_strand,intron_strand_type])
	if exon_strand_type == "intergentic" and intron_strand_type == "intergentic":
		return "intergentic","intergentic","intergentic"
	elif exon_strand_type == "single" and intron_strand_type == "intergentic":
		return "single",exon_strand,"exon"
	elif exon_strand_type == "double" and intron_strand_type == "intergentic":
		return "double",exon_strand,"exon"

	elif exon_strand_type == "single" and intron_strand_type == "single":
		return "single",exon_strand,"exon"
	elif exon_strand_type == "double" and intron_strand_type == "single":
		return "double",exon_strand,"exon"
	elif exon_strand_type == "intergentic" and intron_strand_type == "single":
		return "single",intron_strand,"intron"

	elif exon_strand_type == "single" and intron_strand_type == "double":
		return "single",exon_strand,"exon"
	elif exon_strand_type == "double" and intron_strand_type == "double":
		return "double",exon_strand,"exon"
	elif exon_strand_type == "intergentic" and intron_strand_type == "double":
		return "double",intron_strand,"intron"



def exon_strand_estimate(overlap_exon_str):
	out_exon_strand_list = []
	out_exon_strand_dict = {"+":False,"-":False}
	if overlap_exon_str=="False":
		pass	
	else:
		overlap_exon_list = overlap_exon_str.split("#")
		for i in range(len(overlap_exon_list)):
			sub_overlap_exon = overlap_exon_list[i].split(";")
			sub_overlap_exon_strand = sub_overlap_exon[5]
			out_exon_strand_list.append(sub_overlap_exon_strand)
			out_exon_strand_dict[sub_overlap_exon_strand] = True
	strand,strand_type =  exon_snp_strand_estimate(out_exon_strand_list,out_exon_strand_dict)
	return strand,strand_type


def intron_strand_estimate(overlap_intron_str):
	out_intron_strand_list = []
	out_intron_strand_dict = {"+":False,"-":False}
	if overlap_intron_str=="False":
		pass	
	else:
		overlap_intron_list = overlap_intron_str.split("#")
		for i in range(len(overlap_intron_list)):
			sub_overlap_intron = overlap_intron_list[i].split(";")
			sub_overlap_intron_strand = sub_overlap_intron[5]
			out_intron_strand_list.append(sub_overlap_intron_strand)
			out_intron_strand_dict[sub_overlap_intron_strand] = True
	strand,strand_type = intron_snp_strand_estimate(out_intron_strand_list,out_intron_strand_dict)
	return strand,strand_type
	

def alu_strand_estimate(overlap_alu_str):
	out_alu_strand_list = []
	out_alu_strand_dict = {"+":False,"-":False}
	if overlap_alu_str == "False":
		pass
	else:
		overlap_alu_list = overlap_alu_str.split("#")
		for sub_overlap_alu in overlap_alu_list:
			overlap_alu_infor = sub_overlap_alu.split(";")
			alu_strand = overlap_alu_infor[3]
			out_alu_strand_dict[alu_strand] = True
			out_alu_strand_list.append(alu_strand)
	strand,strand_type = alu_snp_strand_estimate(out_alu_strand_list,out_alu_strand_dict)
	return strand,strand_type

def exon_snp_strand_estimate(out_exon_strand_list,out_exon_strand_dict):
	if out_exon_strand_dict["+"] == True and out_exon_strand_dict["-"] == False:
		return "+","single"
	elif out_exon_strand_dict["+"] == False and out_exon_strand_dict["-"] == True:
		return "-","single"
	elif out_exon_strand_dict["+"] == True and out_exon_strand_dict["-"] == True:
		return out_exon_strand_list[0],"double"
	else:
		return "False","intergentic"

def intron_snp_strand_estimate(out_intron_strand_list,out_intron_strand_dict):
	if out_intron_strand_dict["+"] == True and out_intron_strand_dict["-"] == False:
		return "+","single"
	elif out_intron_strand_dict["+"] == False and out_intron_strand_dict["-"] == True:
		return "-","single"
	elif out_intron_strand_dict["+"] == True and out_intron_strand_dict["-"] == True:
		return out_intron_strand_list[0],"double"
	else:
		return "False","intergentic"

def alu_snp_strand_estimate(out_alu_strand_list,out_alu_strand_dict):
	if out_alu_strand_dict["+"] == True and out_alu_strand_dict["-"] == False:
		return "+","single"
	elif out_alu_strand_dict["+"] == False and out_alu_strand_dict["-"] == True:
		return "-","single"
	elif out_alu_strand_dict["+"] == True and out_alu_strand_dict["-"] == True:
		return out_alu_strand_list[0],"double"
	else:
		return "False","intergentic"

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--annotation_file', required=True, help="alu and gene annotation file")
	parser.add_argument('--result_file', required=True, help="result file")
	parser.add_argument('--sequence_type',required=True, help="sequence type = fr-strand or rf-strand or non-strand")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	annotation_file = args.annotation_file
	result_file = args.result_file
	sequence_type = args.sequence_type
	gene_annotation_file_read(annotation_file,result_file,sequence_type)

if __name__ == "__main__":
	main()



