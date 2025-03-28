import argparse

def cellbarcode_and_umi_select(R1_list):
	"""
	@A00133:429:HLFLLDSX2:1:1101:2157:1000 1:N:0:CGTCCACC+GTCATGAA
	NTTCCATAGGAGTACCTATCTGCCTAAA
	+
	#FFFFFFFFFFFFFFFFFFFFFFFFFFF
	"""
	readname = R1_list[0]
	sequence = R1_list[1]
	strand = R1_list[2]
	quanlity = R1_list[3]
	cellbarcode = sequence[0:16]
	umi = sequence[16:-1]
	return cellbarcode,umi,readname

def readname_make(readname,cellbarcode,umi):
	tmpname1 = readname.split(" 1:")[0]
	tmpname2 = " 1:"+readname.split(" 1:")[1]
	new_readname = ";".join([tmpname1,cellbarcode,umi]) + tmpname2
	return new_readname

def readname_consistent(readname1,readname2):
	"""
	@A00133:429:HLFLLDSX2:1:1101:2157:1000 1:N:0:CGTCCACC+GTCATGAA
	@A00133:429:HLFLLDSX2:1:1101:2157:1000 2:N:0:CGTCCACC+GTCATGAA
	"""
	list1 = readname1.split(" 1:")
	list2 = readname2.split(" 2:")
	for i in range(len(list1)):
		if list1[i] == list2[i]:
			pass
		else:
			return True
	return False

def result_write(d,new_readname,R2_list):
	d.write(new_readname+"\n")
	for i in range(1,4):
		line = R2_list[i]
		d.write(line+"\n")
	
def fastq_fileread(R1_fastq_file,R2_fastq_file,result_fastq_file):
	f1 = open(R1_fastq_file,'r')
	f2 = open(R2_fastq_file,'r')
	d = open(result_fastq_file,'a')
	tmp_list1 = []
	tmp_list2 = []
	for line1 in f1:
		line2 = f2.readline()
		line1 = line1.strip("\n")
		line2 = line2.strip("\n")
		tmp_list1.append(line1)
		tmp_list2.append(line2)
		if len(tmp_list1) == 4:
			if readname_consistent(readname1=tmp_list1[0],readname2=tmp_list2[0]):
				raise Exception("name is not similar")
			cellbarcode,umi,readname = cellbarcode_and_umi_select(R1_list=tmp_list1)
			new_readname = readname_make(readname=readname,cellbarcode=cellbarcode,umi=umi)
			result_write(d=d,new_readname=new_readname,R2_list=tmp_list2)
			################################################################
			tmp_list1 = []
			tmp_list2 = []
		else:
			pass


def make_args():
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-R1', '--R1_fastq', required=True, help="fastq with read")
	parser.add_argument('-R2', '--R2_fastq', required=True, help="fastq with umi and cell barcode")
	parser.add_argument('-O', '--output_fastq', required=True, help="call barcode and umi merged fastq")
	args = parser.parse_args()
	return args

def main():
	import sys
	args = make_args()
	R1_fastq_file = args.R1_fastq
	R2_fastq_file = args.R1_fastq
	result_fastq_file = args.output_fastq
	fastq_fileread(R1_fastq_file=R1_fastq_file,R2_fastq_file=R2_fastq_file,result_fastq_file=result_fastq_file)

if __name__=="__main__":
	main()
