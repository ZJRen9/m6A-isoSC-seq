

def matrix2dict(filename,cellid2celltype_dict):
	f = open(filename,'r')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "transId":
			barcodelist=list_x[1:]
			#print(len(barcodelist))
			continue
		valuelist = [int(x) for x in list_x[1:]]
		snpid = list_x[0]
		celltype2value_dict = cellbarcode_list_to_index(barcodelist=barcodelist,cellid2celltype_dict=cellid2celltype_dict,valuelist=valuelist)
		sum_value_1 = celltype2value_dict["HEK293T"]
		sum_value_2 = celltype2value_dict["Hela"]
		sum_value_3 = celltype2value_dict["HepG2"]
		print("\t".join([snpid,str(sum_value_1),str(sum_value_2),str(sum_value_3)]))

def cellbarcode_list_to_index(barcodelist,cellid2celltype_dict,valuelist):
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
		value_i = valuelist[i]
		##################
		try:
			if value_i >= 1:
				outdict[celltype_i] += 1
			else:
				outdict[celltype_i] += 0
		except:
			if value_i >= 1:
				outdict[celltype_i] = 1
			else:
				outdict[celltype_i] = 0
		##################
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
		barcode = list_x[0].split("-")[0]
		celltype = list_x[1]
		#############
		outdict[barcode] = celltype
		#############
	return outdict

def main():
	import sys
	matrix_filename = sys.argv[1]
	celltype_filename = sys.argv[2]
	cellid2celltype_dict = celltype2dict(filename=celltype_filename)
	matrix2dict(filename=matrix_filename,cellid2celltype_dict=cellid2celltype_dict)

if __name__=="__main__":
	main()


