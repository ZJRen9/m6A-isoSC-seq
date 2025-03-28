import argparse


def merge_transcript_fileread(filename,binSize,slide_lenght,outfilename):
	"""
	ENSG00000243485 chr1    +       29553   30039   29553   30039   3       29553,30266,30975,      30039,30667,31109,      ENSG00000243485     MIR1302-2HG     lncRNA
	ENSG00000237613 chr1    -       34553   35174   34553   35174   3       34553,35244,35720,      35174,35481,36081,      ENSG00000237613     FAM138A lncRNA
	ENSG00000186092 chr1    +       65418   65433   65418   65433   3       65418,65519,69036,      65433,65573,71585,      ENSG00000186092     OR4F5   protein_coding
	"""
	f = open(filename,'r')
	d = open(outfilename,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#####################
		chrom = list_x[1]
		strand = list_x[2]
		geneStart = list_x[3]
		geneEnd = list_x[4]
		cdsStart = list_x[5]
		cdsEnd = list_x[6]
		exonNum = list_x[7]
		exon_start_list = list_x[8].split(",")[:-1]
		exon_end_list = list_x[9].split(",")[:-1]
		exon_start_list = [int(x) for x in exon_start_list]
		exon_end_list = [int(x) for x in exon_end_list]
		geneid = list_x[10]
		genename = list_x[11]
		genetype = list_x[12]
		#####################
		#print(exon_start_list)
		#print(exon_end_list)
		binlist = exonlist_to_list(exon_start_list=exon_start_list,exon_end_list=exon_end_list,binSize=binSize,slide_lenght=slide_lenght)
		#####################
		resultwrite(binlist=binlist,chrom=chrom,strand=strand,geneid=geneid,genename=genename,genetype=genetype,d=d)
		#genelength = translength_count(exon_start_list=exon_start_list,exon_end_list=exon_end_list)
		#####################

def resultwrite(binlist,chrom,strand,geneid,genename,genetype,d):
	for sublist in binlist:
		#print(sublist)
		bin_start = str(sublist[0])
		bin_end = str(sublist[1])
		exon_number = str(sublist[2])
		bin_index = str(sublist[3])
		newname = "_".join([geneid,genename,exon_number,bin_index])
		
		str2write = "\t".join([chrom,bin_start,bin_end,newname,"100",strand]) + "\n"
		d.write(str2write)



def exonlist_to_list(exon_start_list,exon_end_list,binSize,slide_lenght):
	tmpextlength = 0
	binindex = 0
	outlist = []
	for i in range(len(exon_start_list)):
		exon_start_i = exon_start_list[i]
		exon_end_i = exon_end_list[i]
		exon_number = i
		#print(exon_start_i,exon_end_i)
		newoutlist = exonsplit(exon_start_i=exon_start_i,exon_end_i=exon_end_i,exon_number=i,binSize=binSize,slide_lenght=slide_lenght)
		#tmpextlength = newextlength
		#binindex = newbinindex
		outlist+=newoutlist
	return outlist

def exonsplit(exon_start_i,exon_end_i,exon_number,binSize,slide_lenght):
	outlist = []
	tmpstart = exon_start_i
	tmpindex = 0
	while tmpstart + binSize < exon_end_i:
		bin_start = tmpstart
		bin_end = tmpstart + binSize
		outlist.append([bin_start,bin_end,exon_number,tmpindex])
		tmpindex += 1
		tmpstart += slide_lenght
	outlist.append([tmpstart,exon_end_i,exon_number,tmpindex])
	return outlist



def translength_count(exon_start_list,exon_end_list):
	total_lenght = 0
	for i in range(len(exon_start_list)):
		exon_start = exon_start_list[i]
		exon_end = exon_end_list[i]
		total_lenght += exon_end - exon_start
	return total_lenght

def make_args():
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-B', '--input_bin_size', required=True, help="input bin size")
	parser.add_argument('-S', '--input_slide_lenght', required=True, help="input slide lenght")
	parser.add_argument('-I', '--input_merged_genePred', required=True, help="input merged genePred")
	parser.add_argument('-O', '--output_window_filename', required=True, help="output window filename")
	args = parser.parse_args()
	return args

def main():
	import sys
	args = make_args()
	filename = args.input_merged_genePred
	window_filename = args.output_window_filename
	binSize = int(args.input_bin_size)
	slide_lenght = int(args.input_slide_lenght)
	######################
	######################
	merge_transcript_fileread(filename=filename,binSize=binSize,slide_lenght=slide_lenght,outfilename = window_filename)

if __name__=="__main__":
	main()



