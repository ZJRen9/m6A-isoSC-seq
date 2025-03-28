def mtx_file_read(filename):
	"""
	%%MatrixMarket matrix coordinate real general
	% written by sprs
	3069 2700 112625
	1 84 0
	1 180 0
	1 227 0
	"""
	f = open(filename,'r')
	f.readline()
	f.readline()
	f.readline()
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split(" ")
		site_number = list_x[0]
		cell_number = list_x[1]
		value = int(list_x[2])
		try:
			outdict[site_number].append([cell_number,value])
		except:
			outdict[site_number] = [[cell_number,value]]
	return outdict

def cellbarcode_to_dict(cellbarcode_filename):
	"""
	AAACATACAACCAC-1
	AAACATTGAGCTAC-1
	AAACATTGATCAGC-1
	"""
	f = open(cellbarcode_filename,'r')
	cell_number = 0
	cell_dict = {}
	cell_list = []
	for str_x in f:
		cell_number +=1
		cell_barcode = str_x.strip("\n")
		cell_dict[str(cell_number)] = cell_barcode
		cell_list.append(str(cell_number))
	return cell_dict,cell_list

def vcf_to_dict(vcf_filename):
	"""
	21      34988737        A2I_1941        T       C       109.77  PASS    AC=1;       GT:AD:DP:GQ:PL   0/1:225,25:250:99:138,0,7225
	"""
	f = open(vcf_filename,'r')
	snp_num = 0
	snp_dict = {}
	snp_list = []
	for str_x in f:
		if str_x[0] == "#":
			continue
		snp_num += 1
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		snp_id = ";".join([list_x[0],list_x[1],list_x[2],list_x[3],list_x[4],list_x[5]])
		snp_dict[str(snp_num)] = snp_id
		snp_list.append(str(snp_num))
	return snp_dict,snp_list

def level_result_write1(ref_mtx_dict,alt_mtx_dict,cell_dict,cell_list,snp_dict,snp_list,snp_level_filename):
	d = open(snp_level_filename,'a')
	snp_level_str = "\t".join(['chrom','position','snp_name','ref','alt','pvalue','ref_depth','alt_depth','ref_cell_num','alt_cell_num']) + "\n"
	d.write(snp_level_str)
	for snp_num in snp_list:
		snpname = snp_dict[snp_num]
		try:
			ref_value_list = ref_mtx_dict[snp_num]
		except:
			ref_value_list = []
		try:
			alt_value_list = alt_mtx_dict[snp_num]
		except:
			alt_value_list = []
		ref_sum_value = 0
		alt_sum_value = 0
		ref_cell_number = 0
		alt_cell_number = 0
		for subref in ref_value_list:
			cell_num = subref[0]
			cellname = cell_dict[cell_num]
			ref_value = subref[1]
			if ref_value > 0:
				ref_cell_number += 1
			ref_sum_value += ref_value

		for subalt in alt_value_list:
			cell_num = subalt[0]
			cellname = cell_dict[cell_num]
			alt_value = subalt[1]
			if alt_value > 0:
				alt_cell_number +=1
			alt_sum_value += alt_value

		snp_id_list = snpname.split(";")
		chrom = snp_id_list[0]
		site = snp_id_list[1]
		snp_name = snp_id_list[2]
		ref = snp_id_list[3]
		alt = snp_id_list[4]
		pvalue = snp_id_list[5]
		snp_level_str = "\t".join([chrom,site,snp_name,ref,alt,pvalue,str(ref_sum_value),str(alt_sum_value),str(ref_cell_number),str(alt_cell_number)]) + "\n"
		d.write(snp_level_str)

def level_result_write2(ref_mtx_dict,alt_mtx_dict,cell_dict,cell_list,snp_dict,snp_list,snp_level_filename,cell_type_dict,celltypelist):
	d1 = open(snp_level_filename+".Number.counts.txt",'a')
	d2 = open(snp_level_filename+".Number.cells.txt",'a')


	print(celltypelist)

	str_value_ref = "\t".join([x+"_ref" for x in celltypelist])
	str_value_alt = "\t".join([x+"_alt" for x in celltypelist])
	snp_level_str1 = "\t".join(["chrom","site","snp_name","ref","alt","pvalue",str_value_ref,str_value_alt]) + "\n"
	snp_level_str2 = "\t".join(["chrom","site","snp_name","ref","alt","pvalue",str_value_ref,str_value_alt]) + "\n"
	d1.write(snp_level_str1)
	d2.write(snp_level_str2)

	for snp_num in snp_list:
		snp_id = snp_dict[snp_num]
		try:
			ref_value_list = ref_mtx_dict[snp_num]
		except:
			ref_value_list = []
		try:
			alt_value_list = alt_mtx_dict[snp_num]
		except:
			alt_value_list = []
		#####################################################
		ref_valuedict1 = dict(zip(celltypelist,[0 for x in celltypelist]))
		ref_valuedict2 = dict(zip(celltypelist,[0 for x in celltypelist]))
		for subref in ref_value_list:
			cell_num = subref[0]
			cellname = cell_dict[cell_num]
			try:
				cell_type = cell_type_dict[cellname]
			except:
				cell_type = "non_class"
			ref_value = subref[1]
			ref_valuedict1[cell_type] += ref_value
			if ref_value > 0:
				ref_valuedict2[cell_type] +=1

		alt_valuedict1 = dict(zip(celltypelist,[0 for x in celltypelist]))
		alt_valuedict2 = dict(zip(celltypelist,[0 for x in celltypelist]))
		for subalt in alt_value_list:
			cell_num = subalt[0]
			cellname = cell_dict[cell_num]
			try:
				cell_type = cell_type_dict[cellname]
			except:
				cell_type = "non_class"
			alt_value = subalt[1]
			alt_valuedict1[cell_type] += alt_value
			if alt_value > 0:
				alt_valuedict2[cell_type] +=1
		#####################################################
		str_value_ref1 = "\t".join([str(ref_valuedict1[x]) for x in celltypelist])
		str_value_ref2 = "\t".join([str(ref_valuedict2[x]) for x in celltypelist])
		str_value_alt1 = "\t".join([str(alt_valuedict1[x]) for x in celltypelist])
		str_value_alt2 = "\t".join([str(alt_valuedict2[x]) for x in celltypelist])
		#####################################################
		snp_id_list = snp_id.split(";")
		chrom = snp_id_list[0]
		site = snp_id_list[1]
		snp_name = snp_id_list[2]
		ref = snp_id_list[3]
		alt = snp_id_list[4]
		pvalue = snp_id_list[5]
		#####################################################
		snp_level_str1 = "\t".join([chrom,site,snp_name,ref,alt,pvalue,str_value_ref1,str_value_alt1]) + "\n"
		snp_level_str2 = "\t".join([chrom,site,snp_name,ref,alt,pvalue,str_value_ref2,str_value_alt2]) + "\n"
		#####################################################
		d1.write(snp_level_str1)
		d2.write(snp_level_str2)

def cell_type_dict_make(cell_type_filename):
	"""
	AAACCTGGTCAGGACA-1      h409B2_neuroectoderm    409b2   Neuroectoderm
	AAACGGGTCCTACAGA-1      h409B2_neuroectoderm    409b2   Neuroectoderm
	AAACGGGTCTTACCGC-1      h409B2_neuroectoderm    409b2   Neuroectoderm
	"""
	f = open(cell_type_filename,'r')
	outdict = {}
	celltypelist = ["non_class"]
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		cell_barcode = list_x[0]
		cell_type = list_x[1]
		outdict[cell_barcode] = cell_type
		celltypelist.append(cell_type)
	celltypelist = list(set(celltypelist))
	celltypelist.sort()
	return outdict,celltypelist

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--ref_matrix', required=True, help="ref matrix")
	parser.add_argument('--alt_matrix', required=True, help="alt matrix")
	parser.add_argument('--cell_barcode', required=True, help="cell barcode")
	parser.add_argument('--snp_vcf', required=True, help="snp vcf")
	parser.add_argument('--cell_type_splited', required=True, default="False", help="snp vcf")
	parser.add_argument('--snp_level_filename', required=True, help="snp vcf")
	args = parser.parse_args()
	return args


def main():
	args = make_args()
	ref_mtx_filename = args.ref_matrix
	alt_mtx_filename = args.alt_matrix
	cellbarcode_filename = args.cell_barcode
	vcf_filename = args.snp_vcf
	cell_type_splited = args.cell_type_splited
	snp_level_filename = args.snp_level_filename
	ref_mtx_dict = mtx_file_read(ref_mtx_filename)
	alt_mtx_dict = mtx_file_read(alt_mtx_filename)
	cell_dict,cell_list = cellbarcode_to_dict(cellbarcode_filename)
	snp_dict,snp_list = vcf_to_dict(vcf_filename)
	if cell_type_splited == "False":
		level_result_write1(ref_mtx_dict,alt_mtx_dict,cell_dict,cell_list,snp_dict,snp_list,snp_level_filename)
	else:
		cell_type_dict,celltypelist = cell_type_dict_make(cell_type_splited)
		level_result_write2(ref_mtx_dict,alt_mtx_dict,cell_dict,cell_list,snp_dict,snp_list,snp_level_filename,cell_type_dict,celltypelist)


if __name__ == "__main__":
	main()
