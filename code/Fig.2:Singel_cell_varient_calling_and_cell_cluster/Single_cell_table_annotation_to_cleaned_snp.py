

def table_file_read(type_annotation_result,edit_cutoff_result,snp_cutoff_result,min_pvalue,min_alt_depth,min_total_depth,min_mutation_ratio):
	f = open(type_annotation_result,'r')
	d1 = open(edit_cutoff_result,'a')
	d2 = open(snp_cutoff_result,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if list_x[0] == "chrom":
			index_chrom = list_x.index("chrom")
			index_position = list_x.index("position")
			index_name = list_x.index("snp_name")
			index_ref = list_x.index("ref")
			index_alt = list_x.index("alt")
			index_pvalue = list_x.index("pvalue")
			index_ref_depth = list_x.index("ref_depth")
			index_alt_depth = list_x.index("alt_depth")
			index_ref_cell_num = list_x.index("ref_cell_num")
			index_alt_cell_num = list_x.index("alt_cell_num")
			index_exon = list_x.index("exon@intron")
			index_single = list_x.index("single@double")
			str_write = "\t".join(['chrom','position','snp_name','ref','alt','pvalue','ref_depth','alt_depth','ref_cell_num','alt_cell_num','alu_annotation','exon_annotation','intron_annotation','editing_type','exon@intron','single@double']) + "\n"
			d1.write(str_write)
			d2.write(str_write)
			continue
		snp_name = list_x[index_name]
		ref = list_x[index_ref]
		alt = list_x[index_alt]
		pvalue = list_x[index_pvalue]
		ref_depth = list_x[index_ref_depth]
		alt_depth = list_x[index_alt_depth]
		exon = list_x[index_exon]
		single = list_x[index_single]
		if editing_cutoff(snp_name,pvalue,ref,alt,ref_depth,alt_depth,single,min_pvalue,min_alt_depth,min_total_depth,min_mutation_ratio):
			#print(pvalue,ref_depth,alt_depth)
			d1.write(str_x + "\n")
		if snp_cutoff(snp_name,pvalue,ref,alt,ref_depth,alt_depth,single,min_pvalue,min_alt_depth,min_total_depth,min_mutation_ratio):
			d2.write(str_x + "\n")

def editing_cutoff(snp_name,pvalue,ref,alt,ref_depth,alt_depth,single,min_pvalue,min_alt_depth,min_total_depth,min_mutation_ratio):
	#print("######################################################")
	ref_umi_number = int(ref_depth)
	alt_umi_number = int(alt_depth)
	if len(ref)>1 or len(alt)>1:
		return False
	if not snp_name == ".":
		return False
	if ref_umi_number+alt_umi_number <= min_total_depth:
		return False
	if alt_umi_number <= min_alt_depth:
		return False
	if float(pvalue) <= min_pvalue:
		return False
	if alt_umi_number/(ref_umi_number+alt_umi_number)<= min_mutation_ratio:
		return False
	if single == "double":
		return False
	if single == "intergentic":
		return False
	return True

def snp_cutoff(snp_name,pvalue,ref,alt,ref_depth,alt_depth,single,min_pvalue,min_alt_depth,min_total_depth,min_mutation_ratio):
	ref_umi_number = int(ref_depth)
	alt_umi_number = int(alt_depth)
	if len(ref)>1 or len(alt)>1:
		return False
	if snp_name == ".":
		return False
	if ref_umi_number+alt_umi_number < min_total_depth:
		return False
	if alt_umi_number < min_alt_depth:
		return False
	if float(pvalue) < min_pvalue:
		return False
	if alt_umi_number/(ref_umi_number+alt_umi_number)<= min_mutation_ratio:
		return False
	if single == "double":
		return False
	if single == "intergentic":
		return False
	return True

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('--type_annotation_result', required=True, help="type annotation result")
	parser.add_argument('--edit_cutoff_result', required=True, help="edit cutoff result")
	parser.add_argument('--snp_cutoff_result', required=True, help="snp cutoff result")
	parser.add_argument('--min_pvalue', required=True, help="min_pvalue")
	parser.add_argument('--min_alt_depth', required=True, help="min_alt_depth")
	parser.add_argument('--min_total_depth', required=True, help="min_total_depth")
	parser.add_argument('--min_mutation_ratio', required=True, help="min_mutation_ratio")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	type_annotation_result = args.type_annotation_result
	edit_cutoff_result = args.edit_cutoff_result
	snp_cutoff_result = args.snp_cutoff_result
	min_pvalue = float(args.min_pvalue)
	min_alt_depth = float(args.min_alt_depth)
	min_total_depth = float(args.min_total_depth)
	min_mutation_ratio = float(args.min_mutation_ratio)
	#print(min_pvalue,min_alt_depth,min_total_depth)
	table_file_read(type_annotation_result=type_annotation_result,
		            edit_cutoff_result=edit_cutoff_result,
		            snp_cutoff_result=snp_cutoff_result,
		            min_pvalue=min_pvalue,
		            min_alt_depth=min_alt_depth,
		            min_total_depth=min_total_depth,
		            min_mutation_ratio=min_mutation_ratio)

if __name__ == "__main__":
	main()



