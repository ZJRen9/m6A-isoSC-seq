import re
import argparse

def make_args():
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-I1', '--vcf_input', required=True, help="vcf result")
	parser.add_argument('-I2', '--mpileup_input', required=True, help="mpileup result")
	parser.add_argument('-O', '--count_result', required=True, help="merge result")
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

def mpileup_result_read(filename,vcfdict,result_file):
	"""
	#chr1    14590   G
	#49      
	#.$.A......,.A..........A......AA......A..........A      
	#=8>777777B7>7777777777>777777>>777777>7777777774>       
	#]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]       
	#150,133,127,126,126,125,124,124,115,106,101,101,101,100,92,90,86,75,73,73,72,72,66,57,48,48,44,39,38,37,28,29,27,26,26,26,23,23,21,21,22,21,27,21,17,11,11,6,37 
	#A00159:913:HM5JVDSX2:2:1427:20130:21198,A00159:913:HM5JVDSX2:2:1270:22688:32268,A00159:913:HM5JVDSX2:2:2661:16143:9439,A00159:913:HM5JVDSX2:2:2414:21838:11506,A00159:913:HM5JVDSX2:2:2448:3640:25144,A00159:913:HM5JVDSX2:2:2205:1271:30749,A00159:913:HM5JVDSX2:2:1470:14317:28823,A00159:913:HM5JVDSX2:2:2466:17309:3881,A00159:913:HM5JVDSX2:2:1267:2700:16250,A00159:913:HM5JVDSX2:2:2619:16595:22748,A00159:913:HM5JVDSX2:2:2478:18756:33098,A00159:913:HM5JVDSX2:2:2605:32859:29590,A00159:913:HM5JVDSX2:2:1325:2130:33833,A00159:913:HM5JVDSX2:2:1366:31222:1329,A00159:913:HM5JVDSX2:2:2142:2980:5556,A00159:913:HM5JVDSX2:2:2166:21178:4131,A00159:913:HM5JVDSX2:2:2342:22164:6026,A00159:913:HM5JVDSX2:2:1123:18539:27085,A00159:913:HM5JVDSX2:2:1568:4851:22514,A00159:913:HM5JVDSX2:2:2449:2591:26428,A00159:913:HM5JVDSX2:2:1245:17074:20071,A00159:913:HM5JVDSX2:2:1172:9480:17754,A00159:913:HM5JVDSX2:2:1440:29405:5321,A00159:913:HM5JVDSX2:2:2335:24542:36886,A00159:913:HM5JVDSX2:2:1469:12373:8234,A00159:913:HM5JVDSX2:2:1673:1325:26991,A00159:913:HM5JVDSX2:2:2266:32307:18490,A00159:913:HM5JVDSX2:2:1435:26720:25692,A00159:913:HM5JVDSX2:2:2215:7355:8750,A00159:913:HM5JVDSX2:2:1168:5168:32988,A00159:913:HM5JVDSX2:2:1547:13575:23907,A00159:913:HM5JVDSX2:2:2327:18394:14184,A00159:913:HM5JVDSX2:2:2654:14271:12148,A00159:913:HM5JVDSX2:2:1411:25599:33301,A00159:913:HM5JVDSX2:2:1677:9932:28369,A00159:913:HM5JVDSX2:2:2543:4499:14450,A00159:913:HM5JVDSX2:2:2257:12970:2722,A00159:913:HM5JVDSX2:2:1239:4670:21981,A00159:913:HM5JVDSX2:2:1308:7157:36589,A00159:913:HM5JVDSX2:2:1644:27100:3552,A00159:913:HM5JVDSX2:2:2106:4472:24674,A00159:913:HM5JVDSX2:2:2227:26946:29403,A00159:913:HM5JVDSX2:2:2146:23050:23469,A00159:913:HM5JVDSX2:2:2635:6659:31970,A00159:913:HM5JVDSX2:2:2660:25617:32362,A00159:913:HM5JVDSX2:2:1250:1967:12571,A00159:913:HM5JVDSX2:2:2653:5629:8860,A00159:913:HM5JVDSX2:2:2575:26549:32033,A00159:913:HM5JVDSX2:2:1517:18828:13996
	#17     
	#,,aa,,,,,a,,,,,,,        
	#CB00BBBBB0BBBBBBC       
	#]]]]]]]]]]]]]]]]]       
	#131,107,101,77,75,75,51,43,36,36,26,23,23,23,20,20,5
	#A00159:913:HM5JVDSX2:2:1667:9805:34381,A00159:913:HM5JVDSX2:2:2441:15049:30530,A00159:913:HM5JVDSX2:2:1526:4562:27649,A00159:913:HM5JVDSX2:2:2621:7581:8735,A00159:913:HM5JVDSX2:2:1670:26268:29763,A00159:913:HM5JVDSX2:2:2258:12590:22420,A00159:913:HM5JVDSX2:2:1639:26241:17315,A00159:913:HM5JVDSX2:2:1431:15953:25081,A00159:913:HM5JVDSX2:2:1325:2130:33833,A00159:913:HM5JVDSX2:2:1601:8431:4664,A00159:913:HM5JVDSX2:2:2342:22164:6026,A00159:913:HM5JVDSX2:2:1366:31222:1329,A00159:913:HM5JVDSX2:2:1508:29134:9392,A00159:913:HM5JVDSX2:2:1523:23547:25269,A00159:913:HM5JVDSX2:2:1267:2700:16250,A00159:913:HM5JVDSX2:2:2476:23213:11443,A00159:913:HM5JVDSX2:2:2166:21178:4131
	"""
	f = open(filename,'r')
	d = open(result_file,'a')
	d.write("\t".join(['chrom','site','ref','alt','site_name','pvalue','countStr','strandStr','baseStr','quanlity','position','readname']) + "\n")
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
		quanty1 = list_x[5]
		poslist1 = list_x[7].split(",")
		readname1 = list_x[8]
		read2 = list_x[10]
		quanty2 = list_x[11]
		poslist2 = list_x[13].split(",")
		readname2 = list_x[14]
		########################################################
		fdict1,rdict1,baselist1 = read_count(ref,read1,1)
		fdict2,rdict2,baselist2 = read_count(ref,read2,2)
		
		count_str1,strand_str1 = base_count(fdict=fdict1,rdict=rdict1)
		count_str2,strand_str2 = base_count(fdict=fdict2,rdict=rdict2)

		base_str1 = base_str_make(baselist1,ref)
		base_str2 = base_str_make(baselist2,ref)

		base_merge,readname_merge,quanty_merge,position_merge = result_merge_by_readname(base_str1=base_str1,base_str2=base_str2,readname1=readname1,readname2=readname2,
			                                                                             quanty1=quanty1,quanty2=quanty2,poslist1=poslist1,poslist2=poslist2)
		merge_fdict,merge_rdict,merge_baselist = read_count(ref,base_merge,1)
		merge_count_str,merge_strand_str = base_count(fdict=merge_fdict,rdict=merge_rdict)

		str_write = "\t".join([chrom,site,ref,alt,site_name,pvalue,merge_count_str,merge_strand_str,base_merge,quanty_merge,position_merge,readname_merge]) + "\n"
		d.write(str_write)

def result_merge_by_readname(base_str1,base_str2,readname1,readname2,quanty1,quanty2,poslist1,poslist2):
	readname_list1 = readname1.split(",")
	readname_list2 = readname2.split(",")
	if len(readname_list1) != len(base_str1) or  len(readname_list2) != len(base_str2):
		raise Exception("length error")
	for i in range(len(base_str2)):
		base_i = base_str2[i]
		readname_i = readname_list2[i]
		quanty_i = quanty2[i]
		position_i = poslist2[i]
		if readname_i in readname_list1:
			pass
		else:
			if base_i in ['a','t','c','g']:
				base_i = base_i.upper()
			elif base_i in ['A','T','C','G']:
				base_i = base_i.lower()
			else:
				pass
			if readname_i == "*":
				continue
			readname_list1.append(readname_i)
			base_str1 = base_str1 + base_i
			quanty1 = quanty1 + quanty_i
			poslist1.append(position_i)
	outbase_str = ""
	outreadname_list = []
	outquanty = ""
	outposition_list = []
	for i in range(len(base_str1)):
		base_i = base_str1[i]
		readname_i = readname_list1[i]
		quanty_i = quanty1[i]
		position_i = poslist1[i]
		if base_i in ['a','t','c','g','A','T','C','G']:
			outbase_str += base_i
			outreadname_list.append(readname_i)
			outquanty += quanty_i
			outposition_list.append(position_i)
		else:
			pass
	return outbase_str,",".join(outreadname_list),outquanty,",".join(outposition_list)

def base_count(fdict,rdict):
	countlist = []
	reverse_count = 0
	forward_count = 0
	for base in ["A","T","C","G"]:
		sumcounts = fdict[base] + rdict[base]
		countlist.append(str(sumcounts))
		reverse_count += rdict[base]
		forward_count += fdict[base]
	count_str = ";".join(countlist)
	strand_str = str(forward_count) + ";" + str(reverse_count)
	return count_str,strand_str
	

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

def main():
	args = make_args()
	vcf_filename = args.vcf_input
	mpileup_file = args.mpileup_input
	result_file = args.count_result
	vcfdict = vcf2dict(vcf_filename)
	mpileup_result_read(mpileup_file,vcfdict,result_file)

if __name__ == "__main__":
	main()

