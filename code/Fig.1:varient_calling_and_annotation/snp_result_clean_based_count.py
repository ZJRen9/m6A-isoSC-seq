#chr1    155391758       T       C       .       134.77  0;12;7;0        
#tcTcCTCttCtTtCTtCTt     
#A;A;<C<CC=C@C=CA;B;     
#86,82,81,81,72,71,69,69,57,50,45,33,33,32,14,14,10,6,1
#A00836:628:HT7T3DSXY:1:1101:23032:22529;CGGGCATTCGCATCGA;CTGAACACATCG,
#A00836:628:HT7T3DSXY:1:1407:30228:32362;TCCGGGATCGATACCG;TAACGCATTGCT,
#A00836:628:HT7T3DSXY:1:1556:18132:8218;TACGGTATCCCATATC;CAATACCGTGAA,
#A00836:628:HT7T3DSXY:1:1446:32777:21371;AAGTTCGCAACTTCCC;ATAAATACTTGG,

import argparse

def counts_file_read(filename,cleaned_snp_result):
	f = open(filename,'r')
	d = open(cleaned_snp_result,'a')
	str2write = "\t".join(['chrom','position','snp_name','ref','alt','pvalue','ref_count','alt_count','baseSequence','baseQuanlity']) + "\n"
	d.write(str2write)
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			continue
		chrom = list_x[0]
		site = list_x[1]
		ref = list_x[2]
		alt = list_x[3]
		snp_name = list_x[4]
		pvalue = list_x[5]
		ATCG_str = list_x[6]

		strand_str = list_x[7]
		base_str = list_x[8]
		base_qual_str = list_x[9]
		position_str = list_x[10]
		sequence_str = list_x[11]
		#print(position_str)
		if position_str == "":
			continue
		sequence_list = readname_extract(sequence_str=sequence_str)
		position_list = position_extract(position_str=position_str)
		base_qual_list = base_qual_extract(base_qual_str=base_qual_str)
		base_str_list = base_extract(base_str=base_str)	
		ATCG_dict = count2dict(ATCG_str=ATCG_str)
		
		out_base_list,out_qual_list,out_ATCG_dict = sequence_substract_by_quanlity(base_str_list=base_str_list,base_qual_list=base_qual_list)

		#print(out_ATCG_dict)
		#print(ATCG_dict,ref,alt,pvalue,snp_name)
		if cutoff(snp_name=snp_name,pvalue=pvalue,ATCG_dict=out_ATCG_dict,ref=ref,alt=alt):
			#print("\t".join([chrom,site,ref,alt,snp_name,pvalue,ATCG_str]))
			str2write = "\t".join([chrom,site,snp_name,ref,alt,pvalue,str(ATCG_dict[ref]),str(ATCG_dict[alt]),base_str,base_qual_str]) + "\n"
			d.write(str2write)
		else:
			pass

def cutoff(snp_name,pvalue,ATCG_dict,ref,alt):
	#print("######################################################")
	if len(ref)>1 or len(alt)>1:
		return False
	ref_umi_number = int(ATCG_dict[ref.upper()])
	alt_umi_number = int(ATCG_dict[alt.upper()])
	if not snp_name == ".":
		#print("1")
		return False
	if ref == alt:
		return False
	if ref_umi_number+alt_umi_number < 25:
		#print("2")
		return False
	if alt_umi_number/(alt_umi_number+ref_umi_number) <= 0.01:
		#print("3")
		return False
	if alt_umi_number/(alt_umi_number+ref_umi_number) >= 0.99:
		#print("4")
		return False
	if alt_umi_number < 5:
		#print("5")
		return False
	if multi_mutation(ATCG_dict)>2:
		#print("6")
		return False
	if float(pvalue) <20:
		#print("7")
		return False
	#print(ref_umi_number,alt_umi_number,)
	return True


def sequence_substract_by_quanlity(base_str_list,base_qual_list):
	out_base_list = []
	out_qual_list = []
	out_ATCG_dict = {"A":0,"T":0,"C":0,"G":0}
	for i in range(len(base_str_list)):
		base_i = base_str_list[i]
		qual_i = base_qual_list[i]
		if qual_i < 0:
			pass
		else:
			out_base_list.append(base_i)
			out_qual_list.append(qual_i)
	for j in range(len(out_base_list)):
		base_j = out_base_list[j]
		out_ATCG_dict[base_j.upper()] +=1
	return out_base_list,out_qual_list,out_ATCG_dict


def count2dict(ATCG_str):
	ATCG_list = ATCG_str.split(";")
	ATCG_list = [int(x) for x in ATCG_list]
	base_list = ["A","T","C","G"]
	outdict = dict(zip(base_list,ATCG_list))
	return outdict

def readname_extract(sequence_str):
	sequence_list = sequence_str.split(",")
	return sequence_list

def position_extract(position_str):
	position_list = [int(x) for x in position_str.split(",")]
	return position_list

def base_qual_extract(base_qual_str):
	base_qual_list = []
	for i in range(len(base_qual_str)):
		base_qual = base_qual_str[i]
		num_base_qual = ord(base_qual) - 33
		base_qual_list.append(num_base_qual)
	return base_qual_list

def base_extract(base_str):
	base_list = []
	for i in range(len(base_str)):
		base = base_str[i]
		base_list.append(base)
	return base_list

def multi_mutation(ATCG_dict):
	count = 0
	for x in ["A","T","C","G"]:
		if ATCG_dict[x] > 0:
			count +=1
		else:
			pass
	return count

def make_args():
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-I', '--input_snp_result', required=True, help="input snp result")
	parser.add_argument('-O', '--cleaned_snp_result', required=True, help="cleaned snp result")
	args = parser.parse_args()
	return args

if __name__ == "__main__":
	import sys
	args = make_args()
	input_snp_result = args.input_snp_result
	cleaned_snp_result = args.cleaned_snp_result
	counts_file_read(filename=input_snp_result,cleaned_snp_result=cleaned_snp_result)
