

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

def resultwrite(ref_value_dict,alt_value_dict,sitename_list,cellid_list,output_ratio_matrix):
	d = open(output_ratio_matrix,'a')
	str2write = "\t".join(["editingID"]+cellid_list) + "\n"
	d.write(str2write)
	for i in range(len(sitename_list)):
		sitename = sitename_list[i]
		tmplist = [sitename]
		for j in range(len(cellid_list)):
			cellid = cellid_list[j]
			try:
				ref_value = ref_value_dict[(cellid,sitename)]
			except:
				ref_value = 0
			try:
				alt_value = alt_value_dict[(cellid,sitename)]
			except:
				alt_value = 0
			if alt_value+ref_value == 0:
				ratio_value = "0.0"
			else:
				ratio_value = 100*alt_value/(alt_value+ref_value)
				ratio_value = str(ratio_value)
			tmplist.append(ratio_value)
		str2write = "\t".join(tmplist) + "\n"
		d.write(str2write)
		############################

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

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--ref_matrix_filename', required=True, help="ref matrix")
	parser.add_argument('--alt_matrix_filename', required=True, help="alt matrix")
	parser.add_argument('--cell_type_filename', required=True, help="cell type")
	parser.add_argument('--output_ratio_matrix', required=True, help="output_ratio_matrix")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	ref_matrix_filename = args.ref_matrix_filename
	alt_matrix_filename = args.alt_matrix_filename
	cell_type_filename = args.cell_type_filename
	output_ratio_matrix = args.output_ratio_matrix
	#######################
	ref_value_dict,ref_sitename_list,ref_cellid_list = matrix2dict(filename=ref_matrix_filename)
	alt_value_dict,alt_sitename_list,alt_cellid_list = matrix2dict(filename=alt_matrix_filename)
	cellid_list= celltype_fileread(filename=cell_type_filename)
	#######################
	resultwrite(ref_value_dict=ref_value_dict,alt_value_dict=alt_value_dict,sitename_list=ref_sitename_list,cellid_list=cellid_list,output_ratio_matrix=output_ratio_matrix)

if __name__=="__main__":
	main()	


