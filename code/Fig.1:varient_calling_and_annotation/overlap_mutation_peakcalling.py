import argparse
def overlap_count_fileread(filename):
	"""
	chr4    68313356        68313556        ENSG00000083896_YTHDC1_0_297    -       8       635     99
	chr4    68313366        68313566        ENSG00000083896_YTHDC1_0_298    -       9       748     110
	chr4    68313376        68313576        ENSG00000083896_YTHDC1_0_299    -       11      962     150
	chr4    68313386        68313586        ENSG00000083896_YTHDC1_0_300    -       12      1082    160
	"""
	f = open(filename,'r')
	outdict = {}
	outlist = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[0]
		start = list_x[1]
		end = list_x[2]
		windowname = list_x[3]
		strand = list_x[4]
		snpnumber = int(list_x[5])
		refcount = int(list_x[6])
		altcount = int(list_x[7])
		##############################
		geneid = windowname.split("_")[0]
		genename = windowname.split("_")[1]
		exonnumber = windowname.split("_")[2]
		windowid = windowname.split("_")[3]
		##############################
		try:
			outdict[geneid].append([chrom,start,end,windowname,strand,snpnumber,refcount,altcount])
		except:
			outdict[geneid] = [[chrom,start,end,windowname,strand,snpnumber,refcount,altcount]]
			outlist.append(geneid)
	return outdict,outlist

def snpnumber_cutoff(outlist,outdict,min_snpnumber,writefilename):
	d = open(writefilename,'a')
	for geneid in outlist:
		tmplist = outdict[geneid]
		##############################
		chrom = tmplist[0][0]
		windowname =tmplist[0][3]
		strand = tmplist[0][4]
		geneid = windowname.split("_")[0]
		genename = windowname.split("_")[1]
		exonnumber = windowname.split("_")[2]
		##############################
		tmplist = [x for x in tmplist if x[5]>= min_snpnumber]
		win_start_list = [int(x[1]) for x in tmplist]
		win_end_list = [int(x[2]) for x in tmplist]
		win_start_list.sort()
		win_end_list.sort()
		if len(win_start_list) < 1:
			continue
		##############################
		mergelist = bedmerge(win_start_list=win_start_list,win_end_list=win_end_list)
		#mergelist = centerbed(mergelist=mergelist)
		resultwrite(chrom=chrom,strand=strand,mergelist=mergelist,geneid=geneid,genename=genename,exonnumber=exonnumber,d=d)
		##############################

def resultwrite(chrom,strand,mergelist,geneid,genename,exonnumber,d):
	for i in range(len(mergelist)):
		merge_start_i = str(mergelist[i][0])
		merge_end_i = str(mergelist[i][1])
		peakname = "_".join([geneid,genename,exonnumber,str(i)])
		str2write = "\t".join([chrom,merge_start_i,merge_end_i,peakname,strand])
		d.write(str2write)



def bedmerge(win_start_list,win_end_list):
	tmplist = [[win_start_list[0],win_end_list[0]]]
	for i in range(len(win_start_list)):
		win_start_i = win_start_list[i]
		win_end_i = win_end_list[i]
		########################################
		tmp_start = tmplist[-1][0]
		tmp_end = tmplist[-1][1]
		########################################
		overlap_type,overlap_start,overlap_end,overlap_len = bed_overlap_bed(start1=win_start_i,end1=win_end_i,start2=tmp_start,end2=tmp_end)
		if overlap_len > 0:
			new_start = min([tmp_start,win_start_i])
			new_end = max([tmp_end,win_end_i])
			tmplist[-1] = [new_start,new_end]
		else:
			tmplist.append([win_start_i,win_end_i])
	return tmplist

def centerbed(mergelist):
	outlist = []
	for sublist in mergelist:
		tmp_start = sublist[0]
		tmp_end = sublist[1]
		center = (tmp_end - tmp_start)//2 + tmp_start
		center_start = center - 50
		center_end = center + 50
		outlist.append([center_start,center_end])
	return outlist
	

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
		raise Exception("bad overlap")
	return overlap_type,overlap_start,overlap_end,overlap_len

def make_args():
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-I', '--counted_intersect_result', required=True, help="input intersect result")
	parser.add_argument('-M', '--min_C_to_T_mutation_number', required=True, help="input intersect result")
	parser.add_argument('-O', '--counted_peak_result', required=True, help="output windown mutation count")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	import sys
	args = make_args()
	counted_intersect_result = args.counted_intersect_result
	min_C_to_T_mutation_number = int(args.min_C_to_T_mutation_number)
	counted_peak_result = args.counted_peak_result
	outdict,outlist = overlap_count_fileread(filename=counted_intersect_result)
	snpnumber_cutoff(outlist=outlist,outdict=outdict,min_snpnumber=min_C_to_T_mutation_number,d=counted_peak_result)

