
def transcript_annotation_to_dict(filename):
	"""
	transcript_id   chrom   strand  trans_start     trans_end       cds_start       cds_end exon_num  exon_start       exon_end        gene_id gene_name       trans_name      gene_type       stop_code_to_junction_distance     max_exon_length max_internal_exon_length        max_cds_length  max_utr3_length    max_utr5_length bed_cds last_exon_length        canonical_transcript_id
	ENST00000456328 chr1    +       11869   14409   11869   14409   3       11869,12613,13221,      12227,12721,14409, ENSG00000223972 DDX11L1 DDX11L1-202     processed_transcript    1188    1188    1081188    0       0       True    1188    ENST00000450305
	ENST00000450305 chr1    +       12010   13670   12010   13670   6       12010,12179,12613,12975,13221,13453,       12057,12227,12697,13052,13374,13670,    ENSG00000223972 DDX11L1 DDX11L1-201     transcribed_unprocessed_pseudogene 217     217     153     217     0       0       True    217     ENST00000450305
	"""
	f = open(filename,'r')
	out_transcript_to_annotation_dict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "transcript_id":
			continue
		transcript_id = list_x[0]
		out_transcript_to_annotation_dict[transcript_id] = list_x
	f.close()
	return out_transcript_to_annotation_dict

def transcript_annotation_fileread(filename,out_transcript_to_annotation_dict,out_filename):
	f = open(filename,'r')
	d = open(out_filename,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "transcript_id":
			index_transcript_id = list_x.index("transcript_id")
			index_trans_start = list_x.index("trans_start")
			index_trans_end = list_x.index("trans_end")
			index_exon_num = list_x.index("exon_num")
			index_strand = list_x.index("strand")
			index_exon_start = list_x.index("exon_start")
			index_exon_end = list_x.index("exon_end")
			index_cds_start = list_x.index("cds_start")
			index_cds_end = list_x.index("cds_end")
			index_canonical_transcript_id = list_x.index("canonical_transcript_id")
			str2write = "\t".join([str_x,'APA_type','canonical','stopcodon_last_exon_overlap']) + "\t"
			d.write(str2write)
			continue
		transcript_id = list_x[index_transcript_id]
		transcript_start = int(list_x[index_trans_start])
		transcript_end = int(list_x[index_trans_end])
		exon_num = int(list_x[index_exon_num])
		strand = list_x[index_strand]
		exon_start = list_x[index_exon_start]
		exon_end = list_x[index_exon_end]
		cds_start = int(list_x[index_cds_start])
		cds_end = int(list_x[index_cds_end])
		################################
		stopcodon_last_exon_overlap = stop_codon_last_exon_overlap(cds_start=cds_start,cds_end=cds_end,strand=strand,exon_start=exon_start,exon_end=exon_end)
		################################
		canonical_transcript_id = list_x[index_canonical_transcript_id]
		################################
		try:
			tmp_annotation_list = out_transcript_to_annotation_dict[canonical_transcript_id]
			if canonical_transcript_id == transcript_id:
				canonical = "YES"
			else:
				canonical = "NO"
		except:
			APA_type = "bed_canonical_transcript_id"
			canonical = "NO"
			str2write = "\t".join([str_x,APA_type,canonical,stopcodon_last_exon_overlap]) + "\n"
			d.write(str2write)
			continue
		canonical_exon_start = tmp_annotation_list[index_exon_start]
		canonical_exon_end = tmp_annotation_list[index_exon_end]
		canonical_exon_num = int(tmp_annotation_list[index_exon_num])
		################################
		if canonical_exon_num <= 1:
			internal_intron_APA = False
			internal_exon_APA = False
		else:
			internal_intron_APA = internal_intron_APA_calculate(
				canonical_exon_start,
				canonical_exon_end,
				transcript_start,
				transcript_end,
				strand
			)
			internal_exon_APA = internal_exon_APA_calculate(
				canonical_exon_start,
				canonical_exon_end,
				transcript_start,
				transcript_end,
				strand
			)
		if internal_intron_APA == True:
			APA_type = "internal_intron_APA"
		elif internal_exon_APA == True:
			APA_type = "internal_exon_APA"
		else:
			APA_type = "Others"
		str2write = "\t".join([str_x,APA_type,canonical,stopcodon_last_exon_overlap]) + "\t"
		
		d.write(str2write)
		################################

def stop_codon_last_exon_overlap(cds_start,cds_end,strand,exon_start,exon_end):
	exon_start_list = [int(x) for x in exon_start.split(",")[:-1]]
	exon_end_list = [int(x) for x in exon_end.split(",")[:-1]]
	exon_start_list.sort()
	exon_end_list.sort()
	###########################
	if strand == "+":
		last_exon_start = exon_start_list[-1]
		last_exon_end = exon_end_list[-1]
		cds_start_true = min([cds_start,cds_end])
		cds_end_true = max([cds_start,cds_end])
		#print(last_exon_start)
		#print(last_exon_end)
		#print(cds_start)
		#print(cds_end)
		if cds_end_true >= last_exon_start and cds_end_true <= last_exon_end:
			return "True"
		else:
			return "False"
	else:
		last_exon_start = exon_start_list[0]
		last_exon_end = exon_end_list[0]
		cds_start_true = max([cds_start,cds_end])
		cds_end_true = min([cds_start,cds_end])
		#print(last_exon_start)
		#print(last_exon_end)
		#print(cds_start)
		#print(cds_end)
		if cds_end_true >= last_exon_start and cds_end_true <= last_exon_end:
			return "True"
		else:
			return "False"
	###########################



def internal_intron_APA_calculate(canonical_exon_start,canonical_exon_end,transcript_start,transcript_end,strand):
	if strand == "+":
		PAS = max([transcript_end,transcript_start])
	else:
		PAS = min([transcript_end,transcript_start])
	canonical_exon_start_list = [int(x) for x in canonical_exon_start.split(",")[:-1]]
	canonical_exon_end_list = [int(x) for x in canonical_exon_end.split(",")[:-1]]
	canonical_exon_start_list.sort()
	canonical_exon_end_list.sort()
	#########################################
	for i in range(len(canonical_exon_start_list)-1):
		intron_start_i = canonical_exon_end_list[i]
		intron_end_i = canonical_exon_start_list[i+1]
		if PAS >= intron_start_i and PAS < intron_end_i:
			return True
	return False 
	#########################################

def internal_exon_APA_calculate(canonical_exon_start,canonical_exon_end,transcript_start,transcript_end,strand):
	if strand == "+":
		PAS = max([transcript_end,transcript_start])
	else:
		PAS = min([transcript_end,transcript_start])
	canonical_exon_start_list = [int(x) for x in canonical_exon_start.split(",")[:-1]]
	canonical_exon_end_list = [int(x) for x in canonical_exon_end.split(",")[:-1]]
	canonical_exon_start_list.sort()
	canonical_exon_end_list.sort()
	#########################################
	#########################################
	if strand == "+":
		for i in range(len(canonical_exon_start_list)-1):
			exon_start_i = canonical_exon_start_list[i]
			exon_end_i = canonical_exon_end_list[i]
			if PAS >= exon_start_i and PAS < exon_end_i:
				return True
	else:
		for i in range(1,len(canonical_exon_start_list)):
			exon_start_i = canonical_exon_start_list[i]
			exon_end_i = canonical_exon_end_list[i]
			if PAS >= exon_start_i and PAS < exon_end_i:
				return True
	return False

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--GenePred_table', required=True, help="GenePred table")
	parser.add_argument('--Intronic_APA_annotation_transcript', required=True, help="alt matrix")
	args = parser.parse_args()
	return args


def main():
	import sys
	args = make_args()
	GenePred_table = args.GenePred_table
	Intronic_APA_annotation_transcript = args.Intronic_APA_annotation_transcript
	
	out_transcript_to_annotation_dict = transcript_annotation_to_dict(filename=genePred_table)
	transcript_annotation_fileread(filename = genePred_table,out_transcript_to_annotation_dict=out_transcript_to_annotation_dict,out_filename=Intronic_APA_annotation_transcript)

if __name__ == "__main__":
	main()
