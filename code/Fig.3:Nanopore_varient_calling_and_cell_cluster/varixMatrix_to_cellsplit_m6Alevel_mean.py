

def matrix2dict(filename1,filename2,cellid2celltype_dict,output_mean_table):
	f1 = open(filename1,'r')
	f2 = open(filename2,'r')
	d1 = open(output_mean_table,'a')
	str2write = "\t".join(['snpid','HEK293T_mean_level','HeLa_mean_level','HepG2_mean_level','HEK293T_poscell_num','HeLa_poscell_num','HepG2_poscell_num']) + "\n"
	d1.write(str2write)
	#####################################################
	for str_x_1 in f1:
		str_x_2 = f2.readline()
		str_x_1 = str_x_1.strip("\n")
		str_x_2 = str_x_2.strip("\n")
		list_x_1 = str_x_1.split("\t")
		list_x_2 = str_x_2.split("\t")
		if list_x_1[0] == "transId":
			barcodelist_1=list_x_1[1:]
			barcodelist_2=list_x_2[1:]
			if barcodelist_1 == barcodelist_2:
				pass
			else:
				raise exception("error")
			#print(len(barcodelist))
			continue
		valuelist_1 = [int(x) for x in list_x_1[1:]]
		valuelist_2 = [int(x) for x in list_x_2[1:]]
		snpid = list_x_1[0]
		celltype2value_dict = cellbarcode_list_to_index(barcodelist=barcodelist_1,cellid2celltype_dict=cellid2celltype_dict,valuelist_1=valuelist_1,valuelist_2=valuelist_2)
		#print(celltype2value_dict)
		try:
			sum_value_1 = sum(celltype2value_dict["HEK293T"])/len(celltype2value_dict["HEK293T"])
			len_value_1 = len(celltype2value_dict["HEK293T"])
		except:
			sum_value_1 = 0
			len_value_1 = 0
		
		try:	
			sum_value_2 = sum(celltype2value_dict["Hela"])/len(celltype2value_dict["Hela"])
			len_value_2 = len(celltype2value_dict["Hela"])
		except:
			sum_value_2 = 0
			len_value_2 = 0
		
		try:
			sum_value_3 = sum(celltype2value_dict["HepG2"])/len(celltype2value_dict["HepG2"])
			len_value_3 = len(celltype2value_dict["HepG2"])
		except:
			sum_value_3 = 0
			len_value_3 = 0
		#len_value_1 = len(celltype2value_dict["HEK293T"])
		#len_value_2 = len(celltype2value_dict["Hela"])
		#len_value_3 = len(celltype2value_dict["HepG2"])
		
		str2write = "\t".join([snpid,str(sum_value_1),str(sum_value_2),str(sum_value_3),str(len_value_1),str(len_value_2),str(len_value_3)]) + "\n"
		#################
		d1.write(str2write)
		#################

def cellbarcode_list_to_index(barcodelist,cellid2celltype_dict,valuelist_1,valuelist_2):
	"""
	editingID   AAACCCAAGCCTGGAA-1  AAACCCAAGCTCGCAC-1
	"""
	outdict = {}
	for i in range(len(barcodelist)):
		barcode_i = barcodelist[i]
		try:
			celltype_i = cellid2celltype_dict[barcode_i]
		except:
			continue
		value_i_1 = valuelist_1[i]
		value_i_2 = valuelist_2[i]
		value_i = value_i_1 + value_i_2
		if value_i_1 >=1:
			level = value_i_1/value_i
		else:
			continue
		##################
		try:
			if value_i_1 >= 1:
				outdict[celltype_i].append(level)
			else:
				pass
		except:
			if value_i_1 >= 1:
				outdict[celltype_i] = [level]
			else:
				outdict[celltype_i] = []
		##################
	#print(outdict)
	return outdict

def celltype2dict(filename):
	"""
	AAACCCAAGCCTGGAA-1      Hela
	AAACCCAAGCTCGCAC-1      HepG2
	AAACCCACACCCAACG-1      Hela
	AAACCCACAGCGACCT-1      HepG2
	AAACCCAGTGAGATCG-1      HEK293T
	"""
	f = open(filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		#############
		barcode = list_x[0].split("-")[0]
		celltype = list_x[1]
		outdict[barcode] = celltype
		#############
	return outdict

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--ref_matrix_filename', required=True, help="ref matrix")
	parser.add_argument('--alt_matrix_filename', required=True, help="alt matrix")
	parser.add_argument('--cell_type_filename', required=True, help="cell type")
	parser.add_argument('--output_mean_table', required=True, help="output mean table")
	args = parser.parse_args()
	return args


def main():
	import sys
	args = make_args()
	####################################
	matrix_filename_1 = args.ref_matrix_filename
	matrix_filename_2 = args.alt_matrix_filename
	celltype_filename = args.cell_type_filename
	output_mean_table = args.output_mean_table
	####################################
	cellid2celltype_dict = celltype2dict(filename=celltype_filename)
	matrix2dict(filename1=matrix_filename_1,filename2=matrix_filename_2,cellid2celltype_dict=cellid2celltype_dict,output_mean_table=output_mean_table)

if __name__=="__main__":
	main()


