import numpy as np

def canonical_transcript_fileread(filename="/home/ZJRen/Single_cell_m6A/result/transcript/annotation/gencode.v43.canonical_transcript.tsv"):
	"""
	gene_id transcript_id
	ENSG00000290825 ENST00000456328
	ENSG00000223972 ENST00000450305
	ENSG00000227232 ENST00000488147
	"""
	f = open(filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "gene_id":
			continue
		geneid = list_x[0]
		canonical_transid = list_x[1]
		outdict[geneid] = canonical_transid
	return outdict

def intron_retained_isoform_defined(canonical_exon_list,target_exon_list):
	canonical_exon_start_list = [x[0] for x in canonical_exon_list]
	canonical_exon_end_list = [x[1] for x in canonical_exon_list]
	canonical_exon_start_list.sort()
	canonical_exon_end_list.sort()
	for i in range(len(canonical_exon_start_list)-1):
		intron_start_i = canonical_exon_end_list[i]
		intron_end_i = canonical_exon_start_list[i+1]
		for j in range(len(target_exon_list)):
			tmp_target_exon_start_j = target_exon_list[j][0]
			tmp_target_exon_end_j = target_exon_list[j][1]
			target_exon_start_j = min(tmp_target_exon_start_j,tmp_target_exon_end_j)
			target_end_end_j = max(tmp_target_exon_start_j,tmp_target_exon_end_j)
			overlap_type,overlap_start,overlap_end,overlap_len = bed_overlap_bed(start1=intron_start_i,end1=intron_end_i,start2=target_exon_start_j,end2=target_end_end_j)
			if overlap_len >= 200 and overlap_type==4:
				return "intron_retained_isoform"
			else:
				pass
	return "Other"

def coding_sequence_escape_EJC(cds_start_list,cds_end_list):
	outlist = []
	coding_first_exon = abs(cds_start_list[0] - cds_end_list[0]) - 100
	coding_last_exon = abs(cds_start_list[-1] - cds_end_list[-1]) - 100
	if coding_first_exon <= 0:
		coding_first_exon = 0
	if coding_last_exon <= 0:
		coding_last_exon = 0
	######################################
	outlist = [coding_first_exon]
	######################################
	######################################
	for i in range(len(cds_start_list[1:-1])):
		tmp_cds_start = cds_start_list[i]
		tmp_cds_end = cds_end_list[i]
		coding_intern_exon = abs(tmp_cds_end - tmp_cds_start) - 200
		if coding_intern_exon <= 0:
			coding_intern_exon = 0
		outlist.append(coding_intern_exon)
	outlist.append(coding_last_exon)
	return outlist

def coding_internal_exon(exon_start_list,exon_end_list,cds_start,cds_end):
	pass

def annotation_file_to_exon_dict(filename):
	"""
	#trans_id       chrom   strand  trans_start     trans_end       cds_start       cds_end exon_num        exon_start     exon_end        gene_id gene_name       gene_type
	ENST00000456328.2       chr1    +       11869   14409   11869   14409   3       11869,12613,13221,      12227,12721,14409,     ENSG00000223972.5       DDX11L1 lncRNA
	ENST00000450305.2       chr1    +       12010   13670   12010   13670   6       12010,12179,12613,12975,13221,13453,   12057,12227,12697,13052,13374,13670,    ENSG00000223972.5       DDX11L1 transcribed_unprocessed_pseudogene
	"""
	f = open(filename,'r')
	transcript_id_to_exon_dict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		##############
		if list_x[0] == "#trans_id" or list_x[0] == "transcript_id":
			continue
		##############
		#print(str_x)
		##############
		list_x[0] = list_x[0].split(".")[0]
		##############
		transcript_id = list_x[0]
		chrom = list_x[1]
		strand = list_x[2]
		transcript_start = int(list_x[3])
		transcript_end = int(list_x[4])
		cds_start = int(list_x[5])
		cds_end = int(list_x[6])
		exon_number = list_x[7]
		exon_start_list = [int(x) for x in list_x[8].split(",")[:-1]]
		exon_end_list = [int(x) for x in list_x[9].split(",")[:-1]]
		gene_id = list_x[10].split(".")[0]
		gene_name = list_x[11]
		trans_name = list_x[12]
		transcript_type = list_x[13]
		exon_list = [[exon_start_list[i],exon_end_list[i]] for i in range(len(exon_start_list))]
		##########################################
		transcript_id_to_exon_dict[transcript_id] = [exon_list,cds_start,cds_end]
		##########################################
	f.close()
	return transcript_id_to_exon_dict



def fileread(filename,geneid_to_canonical_transid_dict,transcript_id_to_exon_dict,output_file):
	"""
	#trans_id       chrom   strand  trans_start     trans_end       cds_start       cds_end exon_num        exon_start     exon_end        gene_id gene_name       gene_type
	ENST00000456328.2       chr1    +       11869   14409   11869   14409   3       11869,12613,13221,      12227,12721,14409,     ENSG00000223972.5       DDX11L1 lncRNA
	ENST00000450305.2       chr1    +       12010   13670   12010   13670   6       12010,12179,12613,12975,13221,13453,   12057,12227,12697,13052,13374,13670,    ENSG00000223972.5       DDX11L1 transcribed_unprocessed_pseudogene
	"""
	f = open(filename,'r')
	d = open(output_file,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		##############
		if list_x[0] == "#trans_id" or list_x[0] == "transcript_id":
			str2write = "\t".join(list_x+['stop_code_to_junction_distance','canonical_stop_code_to_junction_distance','max_exon_length','max_internal_exon_length','max_cds_length','max_utr3_length','max_utr5_length','bed_cds','last_exon_length',"max_internal_CDS_exon_length","self_retained_intron_transcript_type","sum_length_escape_coding","max_length_escape_coding","canonical_type"]) + "\n"
			d.write(str2write)
			continue
		##############
		#print(str_x)
		##############
		list_x[0] = list_x[0].split(".")[0]
		##############
		transcript_id = list_x[0]
		chrom = list_x[1]
		strand = list_x[2]
		transcript_start = int(list_x[3])
		transcript_end = int(list_x[4])
		cds_start = int(list_x[5])
		cds_end = int(list_x[6])
		exon_number = list_x[7]
		exon_start_list = [int(x) for x in list_x[8].split(",")[:-1]]
		exon_end_list = [int(x) for x in list_x[9].split(",")[:-1]]
		exon_list = [[exon_start_list[i],exon_end_list[i]] for i in range(len(exon_start_list))]
		gene_id = list_x[10].split(".")[0]
		gene_name = list_x[11]
		trans_name = list_x[12]
		transcript_type = list_x[13]
		exon_length_list = exon_length_count(exon_start_list=exon_start_list,exon_end_list=exon_end_list)
		max_exon_length = max(exon_length_list)

		try:
			canonical_transcript_id = geneid_to_canonical_transid_dict[gene_id]
			canonical_exon_list = transcript_id_to_exon_dict[canonical_transcript_id][0]
			canonical_cds_start = transcript_id_to_exon_dict[canonical_transcript_id][1]
			canonical_cds_end = transcript_id_to_exon_dict[canonical_transcript_id][2]
			##############################################################
			self_retained_intron_transcript_type = intron_retained_isoform_defined(canonical_exon_list=canonical_exon_list,target_exon_list=exon_list)
			##############################################################
			if canonical_transcript_id == transcript_id:
				canonical_type = "YES"
			else:
				canonical_type = "NO"
		except:
			canonical_type = "NO"
			self_retained_intron_transcript_type = "bed_canonical"
		###############################
		if transcript_start == cds_start or transcript_end == cds_end:
			bed_cds = "True"
		else:
			bed_cds = "False"
		###############################
		if len(exon_start_list) >= 3:
			max_internal_exon_length = max(exon_length_list[1:-1])
		else:
			max_internal_exon_length = 0
		#############################
		if strand == "+":
			transcript_start_j = min([transcript_start,transcript_end])
			transcript_end_j = max([transcript_start,transcript_end])
			cds_start_j = min([cds_start,cds_end])
			cds_end_j = max([cds_start,cds_end])
			cds_start_list,cds_end_list = target_region_abstract(
			input_start=cds_start_j,
			input_end=cds_end_j,
			exon_start_list=exon_start_list,
			exon_end_list=exon_end_list
			)
			UTR3_start_list,UTR3_end_list = target_region_abstract(
			input_start=transcript_start_j,
			input_end=cds_start_j,
			exon_start_list=exon_start_list,
			exon_end_list=exon_end_list
			)
			UTR5_start_list,UTR5_end_list = target_region_abstract(
			input_start=cds_end_j,
			input_end=transcript_end_j,
			exon_start_list=exon_start_list,
			exon_end_list=exon_end_list
			)
			################################
			stop_code_to_junction_distance = cloest_junction_distance_count(input_pos=cds_end_j,exon_start_list=exon_start_list,exon_end_list=exon_end_list)
			last_exon_length = exon_length_list[-1]

			################################
		else:
			transcript_start_j = max([transcript_start,transcript_end])
			transcript_end_j = min([transcript_start,transcript_end])
			cds_start_j = max([cds_start,cds_end])
			cds_end_j = min([cds_start,cds_end])
			########################################
			cds_start_list,cds_end_list = target_region_abstract(
			input_start=cds_end_j,
			input_end=cds_start_j,
			exon_start_list=exon_start_list,
			exon_end_list=exon_end_list
			)
			UTR3_start_list,UTR3_end_list = target_region_abstract(
			input_start=cds_start_j,
			input_end=transcript_start_j,
			exon_start_list=exon_start_list,
			exon_end_list=exon_end_list
			)
			UTR5_start_list,UTR5_end_list = target_region_abstract(
			input_start=transcript_end_j,
			input_end=cds_end_j,
			exon_start_list=exon_start_list,
			exon_end_list=exon_end_list
			)
			#########################
			stop_code_to_junction_distance = cloest_junction_distance_count(input_pos=cds_end_j,exon_start_list=exon_start_list,exon_end_list=exon_end_list)
			last_exon_length = exon_length_list[0]
		#############################
		try:
			if strand == "+":
				canonical_stop_codon = max(canonical_cds_start,canonical_cds_end)
				canonical_stop_code_to_junction_distance = cloest_junction_distance_count(input_pos=canonical_stop_codon,exon_start_list=exon_start_list,exon_end_list=exon_end_list)
			else:
				canonical_stop_codon = min(canonical_cds_start,canonical_cds_end)
				canonical_stop_code_to_junction_distance = cloest_junction_distance_count(input_pos=canonical_stop_codon,exon_start_list=exon_start_list,exon_end_list=exon_end_list)
		except:
			canonical_stop_code_to_junction_distance = "-1"
		#############################
		#############################
		cds_exon_length_list = exon_length_count(exon_start_list=cds_start_list,exon_end_list=cds_end_list)
		utr3_exon_length_list = exon_length_count(exon_start_list=UTR3_start_list,exon_end_list=UTR3_end_list)
		utr5_exon_length_list = exon_length_count(exon_start_list=UTR5_start_list,exon_end_list=UTR5_end_list)
		if len(cds_exon_length_list) == 0:
			max_cds_length = 0
		else:
			max_cds_length = max(cds_exon_length_list)
		if len(utr3_exon_length_list) == 0:
			max_utr3_length = 0
		else:
			max_utr3_length = max(utr3_exon_length_list)
		if len(utr5_exon_length_list) == 0:
			max_utr5_length = 0
		else:
			max_utr5_length = max(utr5_exon_length_list)
		########################################################################
		if len(exon_start_list) >= 3:
			inter_exon_start_list,inter_exon_end_list = target_region_abstract(input_start=cds_start,input_end=cds_end,exon_start_list=exon_start_list[1:-1],exon_end_list=exon_end_list[1:-1])
			internal_cds_exon_length_list = exon_length_count(exon_start_list=inter_exon_start_list,exon_end_list=inter_exon_end_list)
			if len(internal_cds_exon_length_list) == 0:
				max_internal_CDS_exon_length = 0
			else:
				max_internal_CDS_exon_length = max(internal_cds_exon_length_list)
		else:
			max_internal_CDS_exon_length = 0
		#############################
		if len(cds_start_list) == 1:
			escape_coding_length_list = coding_sequence_escape_EJC(cds_start_list=cds_start_list,cds_end_list=cds_end_list)
			sum_length_escape_coding = escape_coding_length_list[0]
			max_length_escape_coding = max(escape_coding_length_list)
		else:
			escape_coding_length_list = coding_sequence_escape_EJC(cds_start_list=cds_start_list,cds_end_list=cds_end_list)
			sum_length_escape_coding = sum(escape_coding_length_list)
			max_length_escape_coding = max(escape_coding_length_list)

		#############################
		str2write = "\t".join(list_x + [str(stop_code_to_junction_distance),str(canonical_stop_code_to_junction_distance),str(max_exon_length),str(max_internal_exon_length),str(max_cds_length),str(max_utr3_length),str(max_utr5_length),bed_cds,str(last_exon_length),str(max_internal_CDS_exon_length),str(self_retained_intron_transcript_type),str(sum_length_escape_coding),str(max_length_escape_coding),str(canonical_type)]) + "\n"
		#print(str2write)
		d.write(str2write)
		#############################


def target_region_abstract(input_start,input_end,exon_start_list,exon_end_list):
	out_start_list = []
	out_end_list = []
	for i in range(len(exon_start_list)):
		exon_start_i = exon_start_list[i]
		exon_end_i = exon_end_list[i]
		(overlap_type,overlap_start,overlap_end,overlap_len) = bed_overlap_bed(
			start1=input_start,
			end1=input_end,
			start2=exon_start_i,
			end2=exon_end_i
			)
		if overlap_len >= 1:
			out_start_list.append(overlap_start)
			out_end_list.append(overlap_end)
		else:
			pass
	return out_start_list,out_end_list

def cloest_junction_distance_count(input_pos,exon_start_list,exon_end_list):
	out_junction_list = []
	for i in range(len(exon_start_list)):
		exon_start_i = exon_start_list[i]
		exon_end_i = exon_end_list[i]
		###########################################
		if input_pos >= exon_start_i:
			distance_3 = input_pos - exon_start_i
		else:
			distance_3 = exon_start_i - input_pos
		if input_pos >= exon_end_i:
			distance_5 = input_pos - exon_end_i
		else:
			distance_5 = exon_end_i - input_pos
		###########################################
		out_junction_list.append(distance_3)
		out_junction_list.append(distance_5)
		###########################################
	tmp_junction_list = out_junction_list[1:-1]
	if tmp_junction_list == []:
		return np.Inf
	else:
		tmp_junction_list.sort()
		return tmp_junction_list[0]


def exon_length_count(exon_start_list,exon_end_list):
	out_exon_length_list = []
	for i in range(len(exon_start_list)):
		exon_start_i = exon_start_list[i]
		exon_end_i = exon_end_list[i]
		length_i = exon_end_i - exon_start_i
		if length_i < 0:
			raise Exception("ERROR")
		out_exon_length_list.append(length_i)
	return out_exon_length_list
		
def pos_overlap_bed(input_pos,start,end):
	if input_pos >= start and input_pos <= end:
		return True
	else:
		return False

def bed_overlap_bed(start1,end1,start2,end2):
	"""
	input:1000,1200,1100,1300
	output:5,1100,1200,100
	use:overlap_type,overlap_start,overlap_end,overlap_len = bed_overlap_bed(start1,end1,start2,end2)
	"""
	########################################################################          1
	if start1 > start2 and start1 > end2 and end1 > start2 and end1 > end2:
		overlap_type = 1
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	# type1 equal
	elif start1 > start2 and start1 == end2 and end1 > start2 and end1 > end2:
		overlap_type = 1
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	########################################################################          2
	elif start1 > start2 and start1 < end2 and end1 > start2 and end1 > end2:
		overlap_type = 2
		overlap_start = start1
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	########################################################################          3
	elif start1 < start2 and start1 < end2 and end1 > start2 and end1 > end2:
		overlap_type = 3
		overlap_start = start2
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	# type3 equal
	elif start1 == start2 and start1 < end2 and end1 > start2 and end1 > end2:
		overlap_type = 3
		overlap_start = start2
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	elif start1 < start2 and start1 < end2 and end1 > start2 and end1 == end2:
		overlap_type = 3
		overlap_start = start2
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	elif start1 == start2 and start1 < end2 and end1 > start2 and end1 == end2:
		overlap_type = 3
		overlap_start = start2
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	#########################################################################         4
	elif start1 > start2 and start1 < end2 and end1 > start2 and end1 < end2:
		overlap_type = 4
		overlap_start = start1
		overlap_end = end1
		overlap_len = overlap_end - overlap_start
	# type4 equal
	elif start1 == start2 and start1 < end2 and end1 > start2 and end1 < end2:
		overlap_type = 4
		overlap_start = start1
		overlap_end = end1
		overlap_len = overlap_end - overlap_start
	elif start1 > start2 and start1 < end2 and end1 > start2 and end1 == end2:
		overlap_type = 4
		overlap_start = start1
		overlap_end = end1
		overlap_len = overlap_end - overlap_start
	#########################################################################          5
	elif start1 < start2 and start1 < end2 and end1 > start2 and end1 < end2:
		overlap_type = 5
		overlap_start = start2
		overlap_end = end1
		overlap_len = overlap_end - overlap_start
	#########################################################################          6
	elif start1 < start2 and start1 < end2 and end1 < start2 and end1 < end2:
		overlap_type = 6
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	# type6 equal
	elif start1 < start2 and start1 < end2 and end1 == start2 and end1 < end2:
		overlap_type = 6
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	########################################################################           7
	elif start1 == end1 or start2 == end2:
		overlap_type = 7
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	else:
		print(start1,end1,start2,end2)
		raise Exception("bad overlap")
	return overlap_type,overlap_start,overlap_end,overlap_len

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='peaks in transcript relative position')
	parser.add_argument('--canonical_transcript_anno', required=True, help="canonical transcript annotation")
	parser.add_argument('--GenePred_table', required=True, help="GenePred table")
	parser.add_argument('--output_file', required=True, help="result output")
	args = parser.parse_args()
	return args


if __name__ == "__main__":
	import sys
	args = make_args()
	canonical_transcript_anno = args.canonical_transcript_anno
	GenePred_table = args.GenePred_table
	output_file = args.output_file
	#####################################
	#####################################
	geneid_to_canonical_transid_dict = canonical_transcript_fileread(filename=canonical_transcript_anno)
	transcript_id_to_exon_dict = annotation_file_to_exon_dict(filename = GenePred_table)
	fileread(filename=GenePred_table,geneid_to_canonical_transid_dict=geneid_to_canonical_transid_dict,transcript_id_to_exon_dict=transcript_id_to_exon_dict,output_file=output_file)
