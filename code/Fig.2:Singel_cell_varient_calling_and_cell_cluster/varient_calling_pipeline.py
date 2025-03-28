import argparse

def make_args():
	parser = argparse.ArgumentParser(description='prepare bam file for snp calling')
	parser.add_argument('-F', '--Fqfile', required=True, help="Fastq list file")
	parser.add_argument('-P', '--Piared', required=True, help="paired-end data")
	parser.add_argument('-O', '--Outpath', required=True, help="output dir name")
	parser.add_argument("--picard_path", required = False, default = "/home/ZJRen/anaconda3/share/picard-2.22.3-0/picard.jar", help = "picard using path")
	parser.add_argument("--gatk_path", required = False, default = "/home/ZJRen/soft/gatk3.8.0-master/GenomeAnalysisTK.jar", help = "gatk using path")
	parser.add_argument("--star_index", required = False, default = "/home/ZJRen/index.file/mm10_index.file/10X.resouce/STAR_10X_index/", help = "star mapping index")
	parser.add_argument("--star_fasta", required = False, default = "/home/ZJRen/index.file/mm10_index.file/10X.resouce/mm10_10X_genome.fa", help = "fasta using for index make")
	parser.add_argument("--star_gtf", required = False, default = "/home/ZJRen/index.file/mm10_index.file/10X.resouce/mm10_10X_trans.gtf", help = "gtf using for index make")
	parser.add_argument("--vcf_path", required = False, default = "/home/ZJRen/index.file/mm10_index.file/10X.resouce/mm10_10090.vcf")
	args = parser.parse_args()
	return args

def fqfile_to_list(Fqfile,Piared,Outpath):
	f = open(Fqfile,"r")
	outlist1 = []
	outlist2 = []
	#print(Piared)
	if Piared == "True":
		for str_x in f:
			list_x = str_x[:-1].split("\t")
			outlist1.append(list_x[0]+" "+list_x[1])
			outlist2.append(list_x[0].split("_R1.fastq")[0])
	else:
		for str_x in f:
			fastq = str_x[:-1]
			outlist1.append(fastq)
			outlist2.append(fastq.split(".fastq")[0])
	return outlist1,outlist2

def mapping_pipline(fastqlist,bamlist,star_fasta,star_gtf):
	for i in range(len(bamlist)):
		pass1_result = bamlist[i]+"_1pass"
		fastq = fastqlist[i]
		#print("STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 80000000000 --outFilterMatchNminOverLread 0.9 --outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 1 --genomeDir %s --readFilesIn %s --outFileNamePrefix %s" %(star_index,fastq,pass1_result))
		print("STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --outSAMattributes All --limitBAMsortRAM 80000000000 --outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 1 --genomeDir %s --readFilesIn %s --outFileNamePrefix %s" %(star_index,fastq,pass1_result))
		print("sleep 100s")
		pass1_tab = bamlist[i]+"_1passSJ.out.tab"
		pass2_genome = bamlist[i]+"_2pass_genome"
		pass2_result = bamlist[i]+"_2pass"
		print("mkdir %s"%(pass2_genome))
		print("chmod +777 %s"%(pass2_genome))
		print("STAR --runThreadN 10 --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --sjdbFileChrStartEnd %s"%(pass2_genome,star_fasta,star_gtf,pass1_tab))
		print("sleep 100s")
		#print("STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 80000000000 --outFilterMatchNminOverLread 0.9 --outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 1 --genomeDir ./%s/ --readFilesIn %s --outFileNamePrefix %s"%(pass2_genome,fastq,pass2_result))
		print("STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --outSAMattributes All --limitBAMsortRAM 80000000000 --outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 1 --genomeDir ./%s/ --readFilesIn %s --outFileNamePrefix %s"%(pass2_genome,fastq,pass2_result))
		#print("rm -rf ./*pass1*")
		print("wait")
		print("sleep 100s")
		print("echo 'mapping step is over!!!'")
		print("###############################################")

def picard_pipline(bamlist,picard_path,star_fasta):
	for i in range(len(bamlist)):
		pass2_bam = bamlist[i]+"_2passAligned.sortedByCoord.out.bam"
		unqiue_bam = bamlist[i]+"_unique.bam"
		picard_bam = bamlist[i]+"_picard"
		print("samtools view -h -b -q 20 %s >> %s"%(pass2_bam,unqiue_bam))
		print("java -jar %s ReorderSam INPUT=%s OUTPUT=%s_reordered.bam SEQUENCE_DICTIONARY=%s"%(picard_path,unqiue_bam,picard_bam,star_fasta))
		print("java -jar %s AddOrReplaceReadGroups SO=coordinate INPUT=%s_reordered.bam OUTPUT=%s_addrg.bam RGID=%s RGLB=rna RGPL=ILLUMINA RGPU=lane1 RGSM=%s"%(picard_path,picard_bam,picard_bam,picard_bam,picard_bam))
		print("samtools index %s_addrg.bam"%(picard_bam))
		print("java -jar %s MarkDuplicates INPUT=%s_addrg.bam OUTPUT=%s_rmdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=%s_metrics.txt"%(picard_path,picard_bam,picard_bam,picard_bam))
		print("samtools index %s_rmdup.bam"%(picard_bam))
		#print("rm -rf ./*pass2*")
		print("wait")
		print("sleep 300s")
		print("echo 'picard step is over!!!'")
		print("###############################################")

def gatk_pipline(bamlist,gatk_path,star_fasta,vcf_path):
	for i in range(len(bamlist)):
		picard_bam = bamlist[i]+"_picard"
		gatk_bam = bamlist[i]+"_gatk"
		print("/usr/bin/java -Xmx16g -jar %s -T SplitNCigarReads -R %s -I %s_rmdup.bam -o %s_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"%(gatk_path,star_fasta,picard_bam,gatk_bam))
		print("/usr/bin/java -Xmx16g -jar %s -T BaseRecalibrator -R %s -I %s_split.bam -o %s_report.grp -knownSites %s"%(gatk_path,star_fasta,gatk_bam,gatk_bam,vcf_path))
		print("/usr/bin/java -Xmx16g -jar %s -T PrintReads -R %s -I %s_split.bam -BQSR %s_report.grp -o %s_recalibration.bam -U ALLOW_N_CIGAR_READS"%(gatk_path,star_fasta,gatk_bam,gatk_bam,gatk_bam))
		print("wait")
		print("sleep 300s")
		print("echo 'gatk step is over!!!'")
		print("###############################################")

def snp_discover1(bam_list,gatk_path,star_fasta,vcf_path):
	allbamlist = []
	for i in range(len(bam_list)):
		gatk_bam = bam_list[i]+"_gatk"
		reca_bam = gatk_bam + "_recalibration.bam"
		out_vcf = bam_list[i]+"_snp.vcf"
		allbamlist.append(reca_bam)
		allsample = " -I ".join(allbamlist)
		print("java -Xmx32g -jar %s -T UnifiedGenotyper -R %s -I %s  --dbsnp %s -o %s -stand_call_conf 0 --genotyping_mode DISCOVERY -nt 1" % (gatk_path,star_fasta,reca_bam,vcf_path,out_vcf))
	print("sleep 300s")
	print("echo 'gatk step is over!!!'")
	print("wait")
	print("################################################")

def snp_discover2(bam_list,gatk_path,star_fasta,vcf_path):
	allbamlist = []
	for i in range(len(bam_list)):
		gatk_bam = bam_list[i]+"_gatk"
		reca_bam = gatk_bam + "_recalibration.bam"
		allbamlist.append(reca_bam)
	allsample = " -I ".join(allbamlist)
	out_vcf = "out_snp.vcf"
	print("java -Xmx32g -jar %s -T UnifiedGenotyper -R %s -I %s  --dbsnp %s -o %s -stand_call_conf 0 --genotyping_mode DISCOVERY -nt 1" % (gatk_path,star_fasta,allsample,vcf_path,out_vcf))
	print("sleep 300s")
	print("echo 'gatk step is over!!!'")
	print("wait")
	print("################################################")

def samtools_mpileup11(bam_list,star_fasta):
	allbamlist = []
	for i in range(len(bam_list)):
		gatk_bam = bam_list[i]+"_gatk"
		reca_bam = gatk_bam + "_recalibration.bam"
		out_vcf = bam_list[i]+"_snp.vcf"
		out_pipeup = bam_list[i]+"_snp.pileup"
		allbamlist.append(reca_bam)
		print("samtools mpileup -B -d 100000 -f %s -l %s -q 30 -Q 17 -o %s %s" % (star_fasta,out_vcf,out_pipeup,reca_bam))
	print("sleep 300s")
	print("echo 'gatk step is over!!!'")
	print("wait")
	print("################################################")

def samtools_mpileup12(bam_list,star_fasta):
	allbamlist = []
	for i in range(len(bam_list)):
		gatk_bam = bam_list[i]+"_gatk"
		reca_bam = gatk_bam + "_recalibration.bam"
		reca_bam1 = gatk_bam + "_recalibration.first.bam"
		reca_bam2 = gatk_bam + "_recalibration.second.bam"
		print("samtools view -h -b -f 64 %s >> %s"%(reca_bam,reca_bam1))
		print("samtools view -h -b -f 128 %s >> %s"%(reca_bam,reca_bam2))
		out_vcf = bam_list[i]+"_snp.vcf"
		out_pipeup = bam_list[i]+"_snp.pileup"
		allbamlist.append(reca_bam)
		print("samtools mpileup -B -d 100000 -f %s -l %s -q 30 -Q 17 -o %s %s %s" % (star_fasta,out_vcf,out_pipeup,reca_bam1,reca_bam2))
	print("sleep 30s")
	print("echo 'gatk step is over!!!'")
	print("wait")
	print("################################################")

def samtools_mpileup21(bam_list,star_fasta):
	allbamlist = []
	for i in range(len(bam_list)):
		gatk_bam = bam_list[i]+"_gatk"
		reca_bam = gatk_bam + "_recalibration.bam"
		out_vcf = bam_list[i]+"_snp.vcf"
		out_pipeup = bam_list[i]+"_snp.pileup"
		allbamlist.append(reca_bam)
	allsample = " ".join(allbamlist)
	out_vcf = "out_snp.vcf"
	out_pipeup = "out_snp.pileup"
	print("samtools mpileup -B -d 100000 -f %s -l %s -q 30 -Q 17 -o %s %s" % (star_fasta,out_vcf,out_pipeup,allsample))
	print("sleep 30s")
	print("echo 'gatk step is over!!!'")
	print("wait")
	print("################################################")

def samtools_mpileup2(bam_list,star_fasta):
	allbamlist = []
	for i in range(len(bam_list)):
		gatk_bam = bam_list[i]+"_gatk"
		reca_bam = gatk_bam + "_recalibration.bam"
		reca_bam1 = gatk_bam + "_recalibration.first.bam"
		reca_bam2 = gatk_bam + "_recalibration.second.bam"
		print("samtools view -h -b -f 64 %s >> %s"%(reca_bam,reca_bam1))
		print("samtools view -h -b -f 128 %s >> %s"%(reca_bam,reca_bam2))
		out_vcf = bam_list[i]+"_snp.vcf"
		out_pipeup = bam_list[i]+"_snp.pileup"
		allbamlist.append(reca_bam1,reca_bam2)
	allsample = " ".join(allbamlist)
	out_vcf = "out_snp.vcf"
	out_pipeup = "out_snp.pileup"
	print("samtools mpileup -B -d 100000 -f %s -l %s -q 30 -Q 17 -o %s %s" % (star_fasta,out_vcf,out_pipeup,allsample))
	print("sleep 30s")
	print("echo 'gatk step is over!!!'")
	print("wait")
	print("################################################")

if __name__=="__main__":
	import sys
	args = make_args()
	star_index = args.star_index
	star_fasta = args.star_fasta
	star_gtf = args.star_gtf
	picard_path = args.picard_path
	gatk_path = args.gatk_path
	vcf_path = args.vcf_path
	Fqfile = args.Fqfile
	Piared = args.Piared
	Outpath = args.Outpath
	############################################################
	fastqlist,bamlist = fqfile_to_list(Fqfile,Piared,Outpath)
	mapping_pipline(fastqlist,bamlist,star_fasta,star_gtf)
	picard_pipline(bamlist,picard_path,star_fasta)
	gatk_pipline(bamlist,gatk_path,star_fasta,vcf_path)
	snp_discover1(bamlist,gatk_path,star_fasta,vcf_path)
	#snp_discover2(bamlist,gatk_path,star_fasta,vcf_path)
	#samtools_mpileup11(bamlist,star_fasta)
	#samtools_mpileup12(bamlist,star_fasta)

