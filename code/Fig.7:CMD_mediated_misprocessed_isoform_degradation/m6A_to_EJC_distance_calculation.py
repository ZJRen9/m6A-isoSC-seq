import numpy as np

def transcript_to_dict(filename):
	'''
	#trans_id       chrom   strand  trans_start     trans_end       cds_start       cds_end exon_num        exon_start     exon_end        gene_id gene_name       trans_name      gene_type
	ENST00000456328 chr1    +       11869   14409   11869   14409   3       11869,12613,13221,      12227,12721,14409,     ENSG00000223972 DDX11L1 DDX11L1-202     error
	ENST00000450305 chr1    +       12010   13670   12010   13670   6       12010,12179,12613,12975,13221,13453,  12057,12227,12697,13052,13374,13670,     ENSG00000223972 DDX11L1 DDX11L1-201     error
	ENST00000488147 chr1    -       14404   29570   14404   29570   11      14404,15005,15796,16607,16858,17233,17606,17915,18268,24738,29534,     14501,15038,15947,16765,17055,17368,17742,18061,18366,24891,29570,      ENSG00000227232        WASH7P  WASH7P-201      error
	'''
	f = open(filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if str_x[0]=="#":
			continue
		if list_x[0] == "geneLength":
			continue
		#print(list_x)
		transid = list_x[0]
		chrom = list_x[1]
		strand = list_x[2]
		gene_start = int(list_x[3])
		gene_end = int(list_x[4])
		cds_start = int(list_x[5])
		cds_end = int(list_x[6])
		exon_number  = list_x[7]
		exon_start_list = [int(x) for x in list_x[8].split(",")[:-1]]
		exon_end_list = [int(x) for x in list_x[9].split(",")[:-1]]
		geneid = list_x[11]
		#####################################
		#####################################
		outdict[transid] = [gene_start,gene_end,cds_start,cds_end,exon_start_list,exon_end_list,strand]
		#####################################
	return outdict

def genomePosition_to_transcriptPosition(genome_position,exon_list):
	transcript_position = 0
	for exon in exon_list:
		exon_start = exon[0]
		exon_end = exon[1]
		if genome_position >= exon_end and genome_position > exon_start:
			transcript_position += (exon_end - exon_start)
		elif genome_position < exon_end and genome_position >= exon_start:
			transcript_position += (genome_position - exon_start)
		elif genome_position < exon_end and genome_position < exon_start:
			transcript_position += 0
		else:
			raise Exception("genomePosition_to_transcriptPosition error")
	return transcript_position

def transcriptPosition_to_genomePosition(transcript_position,exon_list):
	genome_position = exon_list[0][0]
	tmp_transcript_position = transcript_position
	for i in range(len(exon_list)):
		exon_start = exon_list[i][0]
		exon_end = exon_list[i][1]
		tmp_genome_position = genome_position + tmp_transcript_position
		if tmp_genome_position > exon_end and tmp_genome_position > exon_start:
			genome_position = exon_list[i+1][0]
			tmp_transcript_position = tmp_transcript_position - (exon_end - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position > exon_start:
			genome_position = tmp_genome_position
			tmp_transcript_position = tmp_transcript_position - (tmp_genome_position - exon_start)
		elif tmp_genome_position == exon_end and tmp_genome_position > exon_start:
			if i == (len(exon_list)-1):
				genome_position = tmp_genome_position
				return genome_position
			else:
				genome_position = exon_list[i+1][0]
				tmp_transcript_position = tmp_transcript_position - (exon_end - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position == exon_start:
			genome_position = tmp_genome_position
			tmp_transcript_position = tmp_transcript_position - (tmp_genome_position - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position < exon_start:
			pass
		else:
			raise Exception("transcriptPosition_to_genomePosition error")
	return genome_position

def cloest_exon_junction_distance(chromStart,chromEnd,exon_list):
	#####################################
	for i in range(len(exon_list)):
		exon_start_i = exon_list[i][0]
		exon_end_i = exon_list[i][1]
		##########
		if chromStart >= exon_start_i and chromEnd <= exon_end_i:
			distance_3 = chromStart - exon_start_i
			distance_5 = exon_end_i - chromStart
			if i == 0:
				distance_3 = np.Inf
			else:
				pass
			if i == len(exon_list)-1:
				distance_5 = np.Inf
			else:
				pass
			return min([distance_3,distance_5])
		else:
			pass
		##########
	return np.Inf

def last_exon_judge(exon_list,strand,chromStart):
	if strand == "+":
		last_exon_start = exon_list[-1][0]
		last_exon_end = exon_list[-1][1]
	else:
		last_exon_start = exon_list[0][0]
		last_exon_end = exon_list[0][1]
	if chromStart >= last_exon_start and chromStart <= last_exon_end:
		return "True"
	else:
		return "False"

def peak_file_read(filename,output_file,transcript_id_to_exon_dict):
	'''
	chrom   genepos transcript_id   gene_id depth   modnum  site_ratio      start   end     strand  gene_name     transcript_type  exon_number
	chr1    1319408 ENST00000545578 ENSG00000127054 25      4       0.16    1319408 1319408 -       INTS11  protein_coding 18
	chr1    1319342 ENST00000545578 ENSG00000127054 29      14      0.4827586206896552      1319342 1319342 -     INTS11   protein_coding  18
	chr1    1311788 ENST00000540437 ENSG00000127054 18      1       0.05555555555555555     1311788 1311788 -     INTS11   protein_coding  19
	'''
	f = open(filename,'r')
	d = open(output_file,'a')
	total_trans_p = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0]=="Chr" or list_x[0]=="chrom":
			index_chrom = list_x.index("chrom")
			index_start = list_x.index("snp_start")
			index_transcript_id = list_x.index("transcript_id")
			index_strand = list_x.index("strand")
			##############################################
			str_write = "\t".join([str_x,'m6A_to_EJC_dis','last_exon_overlap','m6A_to_stop_code']) +"\n"
			d.write(str_write)
			continue
		chrom = list_x[index_chrom]
		start = int(list_x[index_start])
		end = int(list_x[1]) + 1
		transcript_id = list_x[index_transcript_id]
		strand = list_x[index_strand]
		##########################
		gene_start,gene_end,cds_start,cds_end,exon_start_list,exon_end_list,transcript_strand = transcript_id_to_exon_dict[transcript_id]
		##########################
		exon_list = [[exon_start_list[i],exon_end_list[i]] for i in range(len(exon_start_list))]
		##########################################
		chromStart = start + (end - start)//2
		chromEnd = chromStart + 1
		trans_chromStart = genomePosition_to_transcriptPosition(chromStart,exon_list)
		trans_chromEnd = genomePosition_to_transcriptPosition(chromEnd,exon_list)
		##########################################
		trans_first_exon_start = genomePosition_to_transcriptPosition(exon_list[0][0],exon_list)
		trans_first_exon_end = genomePosition_to_transcriptPosition(exon_list[0][1],exon_list)
		trans_last_exon_start = genomePosition_to_transcriptPosition(exon_list[-1][0],exon_list)
		trans_last_exon_end = genomePosition_to_transcriptPosition(exon_list[-1][1],exon_list)
		##########################################
		trans_cdsStart = genomePosition_to_transcriptPosition(int(cds_start),exon_list)
		trans_cdsEnd = genomePosition_to_transcriptPosition(int(cds_end),exon_list)
		##########################################
		last_exon_overlap = last_exon_judge(exon_list=exon_list,strand=strand,chromStart=chromStart)
		##########################################
		if strand == "+":
			m6A_to_stop_code = trans_chromStart - trans_cdsEnd
		else:
			m6A_to_stop_code = (trans_cdsStart - trans_chromStart)
		##########################################
		m6A_to_EJC_dis = cloest_exon_junction_distance(chromStart=chromStart,chromEnd=chromEnd,exon_list=exon_list)
		##########################################
		str_write = "\t".join([str_x,str(m6A_to_EJC_dis),str(last_exon_overlap),str(m6A_to_stop_code)]) +"\n"
		#print(str_write)
		d.write(str_write)
	f.close()
	d.close()

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='peaks in transcript relative position')
	parser.add_argument('--peak_format_file', required=True, help="peak format file")
	parser.add_argument('--transcript_file', required=True, help="transcript file")
	parser.add_argument('--output_file', required=True, help="result output")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	peak_format_file = args.peak_format_file
	transcript_file = args.transcript_file
	output_file = args.output_file
	############################
	transcript_id_to_exon_dict = transcript_to_dict(filename=transcript_file)
	####################################
	peak_file_read(filename=peak_format_file,output_file=output_file,transcript_id_to_exon_dict=transcript_id_to_exon_dict)

if __name__=="__main__":
	main()