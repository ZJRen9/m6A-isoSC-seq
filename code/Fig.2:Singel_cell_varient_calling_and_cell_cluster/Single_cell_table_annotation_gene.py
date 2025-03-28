
def genome_table_read(filename,genome_size_file,binsize):
	'''
	ENST00000641515 chr1    +   65418   71585   65564   70008   3   65418,65519,69036,  65433,65573,71585,  0   ENSG00000186092
	'''
	f = open(filename,'r')
	exon_dict = empty_dictmake(genome_size_file,binsize)
	intron_dict = empty_dictmake(genome_size_file,binsize)
	for str_x in f:
		if str_x[0]=="#":
			continue
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		geneid = list_x[0]
		chrom = list_x[1]
		strand = list_x[2]
		#print(strand)
		if strand == "+":
			genestart = int(list_x[3])
			geneend = int(list_x[4])
			cds_start = int(list_x[5])
			cds_end = int(list_x[6])
		if strand == "-":
			genestart = int(list_x[4])
			geneend = int(list_x[3])
			cds_start = int(list_x[6])
			cds_end = int(list_x[5])
		exon_number = int(list_x[7])
		exon_start_list = [int(x) for x in list_x[8].split(",")[:-1]]
		exon_end_list = [int(x) for x in list_x[9].split(",")[:-1]]
		transid = list_x[11]
		if len(exon_start_list) >= 2:
			intron_start_list = exon_end_list[:-1]
			intron_end_list = exon_start_list[1:]
			exon_dict = bed_format_dictmake(exon_dict,genestart,geneend,exon_start_list,exon_end_list,chrom,strand,geneid,transid,binsize,"exon")
			intron_dict = bed_format_dictmake(intron_dict,genestart,geneend,intron_start_list,intron_end_list,chrom,strand,geneid,transid,binsize,"intron")
		else:
			exon_dict = bed_format_dictmake(exon_dict,genestart,geneend,exon_start_list,exon_end_list,chrom,strand,geneid,transid,binsize,"exon")
			#intron_dict = bed_format_dictmake(intron_dict,intron_start_list,intron_end_list,chrom,strand,geneid,transid,binsize,"intron")
	return exon_dict,intron_dict

def bed_format_dictmake(inputdict,genestart,geneend,exon_start_list,exon_end_list,chrom,strand,geneid,transid,binsize,typeid):
	for i in range(len(exon_start_list)):
		exon_start = exon_start_list[i]
		exon_end = exon_end_list[i]
		exon_start_bin = exon_start//binsize
		exon_end_bin = exon_end//binsize
		for j in range(exon_start_bin,exon_end_bin+1):
			inputdict[(chrom,j)].append([chrom,genestart,geneend,exon_start,exon_end,strand,geneid,transid,strand,typeid])
			#print(chrom,j)
	return inputdict

def empty_dictmake(genome_size_file,binsize):
	f = open(genome_size_file,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[0]
		genome_size = int(list_x[1])
		genome_bin = genome_size//binsize
		for i in range(genome_bin+1):
			outdict[(chrom,i)] = []
	return outdict

def snp_type_fileread(filename,exon_dict,intron_dict,binsize,sequence_type,gene_annotation_file):
	'''
	chr1    1316800 T       C       .       136.77  0;21;15;0
	'''
	f = open(filename,'r')
	d = open(gene_annotation_file,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			chrom_index = list_x.index("chrom")
			position_index = list_x.index("position")
			d.write(str_x+"\t"+"exon_annotation"+"\t"+"intron_annotation"+"\n")
			continue
		chrom = list_x[chrom_index]
		site = int(list_x[position_index])
		site_bin = site//binsize
		exon_list = exon_dict[(chrom,site_bin)]
		intron_list = intron_dict[(chrom,site_bin)]
		overlap_exon_list = site_overlap_bed(site,exon_list,sequence_type)
		overlap_intron_list = site_overlap_bed(site,intron_list,sequence_type)
		#print(overlap_exon_list)
		#print(overlap_intron_list)
		if len(overlap_exon_list) == 0:
			overlap_exon_strand = "False"
		else:
			overlap_exon_list.sort(key=lambda x:x[-1])
			same_exon_list = same_geneid_dropout(overlap_exon_list)
			same_strand_exon_list = same_strand_dropout(same_exon_list)
			overlap_exon_strand = result_extract(same_strand_exon_list)

		if len(overlap_intron_list) == 0:
			overlap_intron_strand = "False"
		else:
			overlap_intron_list.sort(key=lambda x:x[-1])
			same_intron_list = same_geneid_dropout(overlap_intron_list)
			same_strand_intron_list = same_strand_dropout(same_intron_list)
			overlap_intron_strand = result_extract(same_strand_intron_list)
		result_str = str_x+"\t"+overlap_exon_strand+"\t"+overlap_intron_strand+"\n"
		d.write(result_str)


def site_overlap_bed(site,bed_list,sequence_type):
	outlist = []
	for bed in bed_list:
		chrom,genestart,geneend,exon_start,exon_end,strand,geneid,transid,strand,typeid = bed
		if sequence_type == "5_type":
			absdis = abs(site - genestart)
		elif sequence_type == "3_type":
			absdis = abs(site - geneend)
		else:
			pass
		if site < exon_start and site < exon_end:
			pass
		elif site >= exon_start and site < exon_end:
			outlist.append([chrom,genestart,geneend,exon_start,exon_end,strand,geneid,transid,strand,typeid,absdis])
		elif site > exon_start and site >= exon_end:
			pass
		else:
			pass
	return outlist

def same_geneid_dropout(overlap_exon_list):
	outlist = []
	geneid_list = []
	if len(overlap_exon_list) == 0:
		return overlap_exon_list
	else:
		for bed in overlap_exon_list:
			chrom,genestart,geneend,exon_start,exon_end,strand,geneid,transid,strand,typeid,absdis = bed
			if geneid in geneid_list:
				pass
			else:
				outlist.append(bed)
				geneid_list.append(geneid)
	return outlist

def same_strand_dropout(same_exon_list):
	outlist = []
	strand_list = []
	if len(same_exon_list) == 0:
		return overlap_exon_list
	else:
		for bed in same_exon_list:
			chrom,genestart,geneend,exon_start,exon_end,strand,geneid,transid,strand,typeid,absdis = bed
			if strand in strand_list:
				pass
			else:
				outlist.append(bed)
				strand_list.append(strand)
	return outlist

def result_extract(same_strand_exon_list):
	result_list = []
	for i in range(len(same_strand_exon_list)):
		a = ";".join([str(x) for x in same_strand_exon_list[i][:7]])
		result_list.append(a)
	result_str = "#".join(result_list)
	return result_str

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-S', '--snp_file', required=True, help="snp format result")
	parser.add_argument('-T', '--GenePred', required=True, help="GenePred table result")
	parser.add_argument('-G', '--genome_size_file', required=True, help="3_type or 5_type")
	parser.add_argument('-O', '--output_file', required=True, help="output file name")
	parser.add_argument('--sequence_type', required=True, help="3_type or 5_type")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	snp_filename = args.snp_file
	genetable = args.GenePred
	genome_size_filename = args.genome_size_file
	sequence_type = args.sequence_type
	gene_annotation_file = args.output_file
	binsize = 100000
	exon_dict,intron_dict = genome_table_read(genetable,genome_size_filename,binsize)
	snp_type_fileread(snp_filename,exon_dict,intron_dict,binsize,sequence_type,gene_annotation_file)

if __name__=="__main__":
	main()