import re
import argparse

def make_args():
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-I1', '--vcf_input', required=True, help="vcf result")
	parser.add_argument('-I2', '--mpileup_input', required=True, help="mpileup result")
	args = parser.parse_args()
	return args

def read_count(ref,read,strand_type):
	baselist = []
	typelist = []
	i = 0
	while i < len(read):
		base = read[i]
		#print(i,base)
		if base == "*":
			baselist.append(base)
			typelist.append("amabiguous")
			i+=1
			# this is amabiguous
		elif base == "$":
			baselist = baselist[:len(baselist)-1]
			typelist = typelist[:len(typelist)-1]
			baselist.append(read[i-1:i+1])
			typelist.append("stop")
			i+=1
			# this is a stop for one read
		elif base == "^":
			baselist.append(read[i:i+3])
			typelist.append("start")
			i+=3
			# this is a start for one read
		elif base == "+":
			baselist = baselist[:len(baselist)-1]
			typelist = typelist[:len(typelist)-1]
			num=re.match(r"(\d+)",read[i+1:]).group(0)
			baselist.append(read[i-1:i+(len(num)+int(num))])
			typelist.append("indel")
			i += (len(num)+int(num)+1)
			# this is a indel
		elif base == "-":
			baselist = baselist[:len(baselist)-1]
			typelist = typelist[:len(typelist)-1]
			num=re.match(r"(\d+)",read[i+1:]).group(0)
			baselist.append(read[i-1:i+(len(num)+int(num))])
			typelist.append("delete")
			i += (len(num)+int(num)+1)
			# this is a delete
		elif base == ",":
			baselist.append(",")
			typelist.append("ref reverse")
			i+=1
			# ref reverse
		elif base == "a":
			baselist.append("a")
			typelist.append("alt reverse")
			i+=1
			# alt reverse
		elif base == "t":
			baselist.append("t")
			typelist.append("alt reverse")
			i+=1
			# alt reverse
		elif base == "c":
			baselist.append("c")
			typelist.append("alt reverse")
			i+=1
			# alt reverse
		elif base == "g":
			baselist.append("g")
			typelist.append("alt reverse")
			i+=1
			# alt reverse
		elif base == ".":
			baselist.append(".")
			typelist.append("ref forward")
			i+=1
			# ref foward
		elif base == "A":
			baselist.append("A")
			typelist.append("alt forward")
			i+=1
			# alt forward
		elif base == "T":
			baselist.append("T")
			typelist.append("alt forward")
			i+=1
			# alt forward
		elif base == "C":
			baselist.append("C")
			typelist.append("alt forward")
			i+=1
			# alt forward
		elif base == "G":
			baselist.append("G")
			typelist.append("alt forward")
			i+=1
			# alt forward
		else:
			pass
	if read == "*":
		fdict = {"A":0,"T":0,"C":0,"G":0}
		rdict = {"A":0,"T":0,"C":0,"G":0}
		outbase = ["*"]
	else:
		fdict,rdict,outbase = mutil_result_count(baselist,typelist,ref,strand_type)
	return fdict,rdict,outbase

def mutil_result_count(baselist,typelist,ref,strand_type):
	fdict = {"A":0,"T":0,"C":0,"G":0}
	rdict = {"A":0,"T":0,"C":0,"G":0}
	outbase = []
	for i in range(len(baselist)):
		base_i = baselist[i]
		type_i = typelist[i]
		if type_i == "amabiguous":
			pass
		elif type_i == "stop":
			base_i = base_i[0]
			fdcit,rdict = dict_add(fdict,rdict,base_i,ref,strand_type)

		elif type_i == "start":
			base_i = base_i[-1]
			fdcit,rdict = dict_add(fdict,rdict,base_i,ref,strand_type)

		elif type_i == "indel":
			base_i = base_i[0]
			fdcit,rdict = dict_add(fdict,rdict,base_i,ref,strand_type)

		elif type_i == "delete":
			base_i = base_i[0]
			fdcit,rdict = dict_add(fdict,rdict,base_i,ref,strand_type)

		elif type_i == "ref reverse":
			fdcit,rdict = dict_add(fdict,rdict,base_i,ref,strand_type)

		elif type_i == "ref forward":
			fdcit,rdict = dict_add(fdict,rdict,base_i,ref,strand_type)

		elif type_i == "alt reverse":
			fdcit,rdict = dict_add(fdict,rdict,base_i,ref,strand_type)

		elif type_i == "alt forward":
			fdcit,rdict = dict_add(fdict,rdict,base_i,ref,strand_type)
		else:
			pass
		outbase.append(base_i)
	return fdict,rdict,outbase

def dict_add(fdict,rdict,base_i,ref,strand_type):
	if strand_type == 1:
		if base_i == ".":
			fdict[ref]+=1
		elif base_i == "A":
			fdict[base_i] +=1
		elif base_i == "T":
			fdict[base_i] +=1
		elif base_i == "C":
			fdict[base_i] +=1
		elif base_i == "G":
			fdict[base_i] +=1
		##########################
		elif base_i == ",":
			rdict[ref]+=1
		elif base_i == "a":
			rdict[base_i.upper()]+=1
		elif base_i == "t":
			rdict[base_i.upper()]+=1
		elif base_i == "c":
			rdict[base_i.upper()]+=1
		elif base_i == "g":
			rdict[base_i.upper()]+=1
	elif strand_type == 2:
		if base_i == ".":
			rdict[ref]+=1
		elif base_i == "A":
			rdict[base_i] +=1
		elif base_i == "T":
			rdict[base_i] +=1
		elif base_i == "C":
			rdict[base_i] +=1
		elif base_i == "G":
			rdict[base_i] +=1
		##########################
		elif base_i == ",":
			fdict[ref]+=1
		elif base_i == "a":
			fdict[base_i.upper()]+=1
		elif base_i == "t":
			fdict[base_i.upper()]+=1
		elif base_i == "c":
			fdict[base_i.upper()]+=1
		elif base_i == "g":
			fdict[base_i.upper()]+=1
	else:
		pass
	return fdict,rdict


def mpileup_result_read(filename,vcfdict):
	"""
	chr1    3319257 C   22  .....T.........T...T..  CA@DDACCEECBDDC@CCB==C  ]]]]]]]]]]]]]]]]]]]]]]  
	91,88,85,68,68,67,51,41,37,37,35,34,24,19,16,12,11,11,10,10,6,3 
	D00353:230:HVV7KBCXY:2:1205:5987:24609;TTTATGCCAAGCTGTT;AGGCATCGGT,D00353:230:HVV7KBCXY:2:2201:13314:51734;CGATTGATCAGTTGAC;ATTTAGATGG
	"""
	f = open(filename,'r')
	for str_x in f:
		str_x = str_x.rstrip('\n\r')
		list_x = str_x.split("\t")
		chrom = list_x[0]
		site = list_x[1]
		ref = list_x[2].upper()
		try:
			alt,site_name,pvalue = vcfdict[(chrom,site,ref)]
		except:
			continue
		read1 = list_x[4]
		read2 = list_x[5]
		poslist = list_x[7]
		readname = list_x[8]
		fdict1,rdict1,baselist = read_count(ref,read1,1)
		countlist = []
		reverse_count = 0
		forward_count = 0
		for base in ["A","T","C","G"]:
			sumcounts = fdict1[base] + rdict1[base]
			countlist.append(str(sumcounts))
			reverse_count += rdict1[base]
			forward_count += fdict1[base]
		count_str = ";".join(countlist)
		strand_str = str(forward_count) + ";" + str(reverse_count)
		base_str = base_str_make(baselist,ref)
		outbase,outquanlity,outpos,outreadname = amabiguous_base_remove(base_str=base_str,quanlity_str=read2,pos_str=poslist,readname_str=readname)

		print("\t".join([chrom,site,ref,alt,site_name,pvalue,count_str,strand_str,outbase,outquanlity,outpos,outreadname]))

def base_str_make(baselist,ref):
	outlist = []
	for x in baselist:
		if x == "*":
			outlist.append("*")
		elif x == ",":
			outlist.append(ref.lower())
		elif x == "a":
			outlist.append("a")
		elif x == "t":
			outlist.append("t")
		elif x == "c":
			outlist.append("c")
		elif x == "g":
			outlist.append("g")
		elif x == ".":
			outlist.append(ref.upper())
		elif x == "A":
			outlist.append("A")
		elif x == "T":
			outlist.append("T")
		elif x == "C":
			outlist.append("C")
		elif x == "G":
			outlist.append("G")
	outstr = "".join(outlist)
	return outstr

def amabiguous_base_remove(base_str,quanlity_str,pos_str,readname_str):
	poslist = pos_str.split(",")
	readname_list = readname_str.split(",")
	outbase_list = []
	outpos_list = []
	outquanlity_list = []
	outreadname_list = []
	for i in range(len(base_str)):
		base_i = base_str[i]
		quanlity_i = quanlity_str[i]
		pos_i = poslist[i]
		readname_i = readname_list[i]
		if base_i == "*":
			pass
		else:
			outbase_list.append(base_i)
			outquanlity_list.append(quanlity_i)
			outpos_list.append(pos_i)
			outreadname_list.append(readname_i)
	outbase = "".join(outbase_list)
	outpos = ",".join(outpos_list)
	outquanlity = "".join(outquanlity_list)
	outreadname = ",".join(outreadname_list)
	return outbase,outquanlity,outpos,outreadname
	

def vcf2dict(vcf_filename):
	"""
	chr12   67059381        rs1132252178    T       G       80.28   .       AC=2;AF=1.00;AN=2;DB;
	"""
	f = open(vcf_filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.rstrip('\n\r')
		if str_x[0]=="#":
			continue
		list_x = str_x.split("\t")
		chrom = list_x[0]
		site = list_x[1]
		site_name = list_x[2]
		ref = list_x[3]
		alt = list_x[4]
		pvalue = list_x[5]
		outdict[(chrom,site,ref)] = [alt,site_name,pvalue]
	return outdict

if __name__=="__main__":
	args = make_args()
	vcf_filename = args.vcf_input
	mpileup_file = args.mpileup_input
	vcfdict = vcf2dict(vcf_filename)
	mpileup_result_read(mpileup_file,vcfdict)
	#############################

