


def matrix2dict(filename):
	f = open(filename,'r')
	value_dict = {}
	sitename_list = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "editingID":
			cellid_list = list_x[1:]
			continue
		sitename = list_x[0]
		sitename_list.append(sitename)
		for i in range(len(cellid_list)):
			cellid = cellid_list[i]
			value = int(list_x[i+1])
			if value == 0:
				pass
			else:
				value_dict[(cellid,sitename)] = value
	return value_dict,sitename_list,cellid_list

def celltype_fileread(filename):
	"""
	cellid  cellType
	AAACCCAAGCCTGGAA-1      Hela
	AAACCCAAGCTCGCAC-1      HepG2
	AAACCCACACCCAACG-1      Hela
	AAACCCACAGCGACCT-1      HepG2
	"""
	f = open(filename,'r')
	cellid_list = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "cellid":
			continue
		cellid_list.append(list_x[0])
	return cellid_list

def snpedit_to_dict(filename):
	"""
	chr10;254368;.;C;T;100
	chr10;812182;.;G;A;100
	chr1;853896;853897;GAACC;f      chr1_853904_853905_f
	chr1;853896;853897;GAACC;f      chr1_853922_853923_f
	chr1;965488;965489;AGACC;f      chr1_965490_965491_f
	chr1;1217109;1217110;AAACT;r    chr1_1217082_1217083_r
	"""
	f = open(filename,'r')
	outdict = {}
	outlist = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		miCLIPid = list_x[0]
		editid = list_x[1]
		########################
		sublist = editid.split("_")
		if sublist[3] == "f":
			newedit_name = ";".join([sublist[0],sublist[1],".;C;T;100"])
		else:
			newedit_name = ";".join([sublist[0],sublist[1],".;G;A;100"])

		try:
			outdict[miCLIPid].append(newedit_name)
		except:
			outdict[miCLIPid] = [newedit_name]
			outlist.append(miCLIPid)
		########################
		########################
	return outdict,outlist

def matrix_merge(ref_value_dict,alt_value_dict,miCLIPid_list,miCLIPid2editid_dict,cellid_list,ref_result_filename,alt_result_filename):
	d1 = open(ref_result_filename,'a')
	d2 = open(alt_result_filename,'a')
	ref_result_str = "\t".join(['editingID']+cellid_list) + "\n"
	alt_result_str = "\t".join(['editingID']+cellid_list) + "\n"
	d1.write(ref_result_str)
	d2.write(alt_result_str)
	#######################################
	for miCLIPid in miCLIPid_list:
		editid_list = miCLIPid2editid_dict[miCLIPid]
		ref_result_list = [miCLIPid]
		alt_result_list = [miCLIPid]
		for cellid in cellid_list:
			sum_ref_value = 0
			sum_alt_value = 0
			for editid in editid_list:
				try:
					tmp_ref_value = ref_value_dict[(cellid,editid)]
				except:
					tmp_ref_value = 0
				try:
					tmp_alt_value = alt_value_dict[(cellid,editid)]
				except:
					tmp_alt_value = 0
				sum_ref_value += tmp_ref_value
				sum_alt_value += tmp_alt_value
			ref_result_list.append(str(sum_ref_value))
			alt_result_list.append(str(sum_alt_value))
		ref_result_str = "\t".join(ref_result_list) + "\n"
		alt_result_str = "\t".join(alt_result_list) + "\n"
		d1.write(ref_result_str)
		d2.write(alt_result_str)

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--ref_matrix_filename', required=True, help="ref matrix")
	parser.add_argument('--alt_matrix_filename', required=True, help="alt matrix")
	parser.add_argument('--cell_type_filename', required=True, help="cell type")
	parser.add_argument('--ref_matrix_resultname', required=True, help="ref matrix")
	parser.add_argument('--alt_matrix_resultname', required=True, help="alt matrix")
	parser.add_argument('--miCLIPid2editid_filename', required=True, help="alt matrix")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	ref_matrix_filename = args.ref_matrix_filename
	alt_matrix_filename = args.alt_matrix_filename
	ref_matrix_resultname = args.ref_matrix_resultname
	alt_matrix_resultname = args.alt_matrix_resultname
	cell_type_filename = args.cell_type_filename
	miCLIPid2editid_filename = args.miCLIPid2editid_filename
	#######################
	ref_value_dict,ref_sitename_list,ref_cellid_list = matrix2dict(filename=ref_matrix_filename)
	print("step1")
	alt_value_dict,alt_sitename_list,alt_cellid_list = matrix2dict(filename=alt_matrix_filename)
	print("step2")
	cellid_list= celltype_fileread(filename=cell_type_filename)
	print("step3")
	miCLIPid2editid_dict,miCLIPid_list = snpedit_to_dict(filename=miCLIPid2editid_filename)
	print("step4")
	#######################
	matrix_merge(ref_value_dict=ref_value_dict,
				alt_value_dict=alt_value_dict,
				miCLIPid_list=miCLIPid_list,
				miCLIPid2editid_dict=miCLIPid2editid_dict,
				cellid_list=cellid_list,
				ref_result_filename=ref_matrix_resultname,
				alt_result_filename=alt_matrix_resultname)
	print("step5")
	#resultwrite(ref_value_dict=ref_value_dict,alt_value_dict=alt_value_dict,sitename_list=ref_sitename_list,cellid_list=cellid_list)

if __name__=="__main__":
	main()
