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
		site_number = int(list_x[0])
		cell_number = int(list_x[1])
		value = int(list_x[2])
		outdict[(cell_number,site_number)] = value
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
		cell_list.append(str(cell_barcode))
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
		snp_list.append(str(snp_id))
	return snp_dict,snp_list


def resultwrite(mtxdict,cell_list,cell_dict,snp_list,snp_dict,output_filename):
	d = open(output_filename,'a')
	###################################
	str2write = "\t".join(["editingID"]+cell_list) + "\n"
	d.write(str2write)
	###################################
	for snpNumber in range(1,len(snp_list)+1):
		snpid = snp_dict[str(snpNumber)]
		tmplist = [snpid]
		for cellNumber in range(1,len(cell_list)+1):
			try:
				value = mtxdict[(cellNumber,snpNumber)]
			except:
				value = 0
			tmplist.append(str(value))
		str2write = "\t".join(tmplist) + "\n"
		d.write(str2write)

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--matrix_filename', required=True, help="matrix")
	parser.add_argument('--barcode_filename', required=True, help="barcode")
	parser.add_argument('--vcf_filename', required=True, help="vcf")
	parser.add_argument('--output_matrix', required=True, help="output ref or alt matrix filename")
	args = parser.parse_args()
	return args


def main():
	args = make_args()
	matrix_file = args.matrix_filename
	barcode_file = args.barcode_filename
	vcf_file = args.vcf_filename
	output_filename = args.output_matrix
	mtxdict = mtx_file_read(filename=matrix_file)
	cell_dict,cell_list = cellbarcode_to_dict(cellbarcode_filename=barcode_file)
	snp_dict,snp_list = vcf_to_dict(vcf_filename=vcf_file)
	resultwrite(mtxdict=mtxdict,cell_list=cell_list,cell_dict=cell_dict,snp_list=snp_list,snp_dict=snp_dict,output_filename=output_filename)

if __name__=="__main__":
	main()
