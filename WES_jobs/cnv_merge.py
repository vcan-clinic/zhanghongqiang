#-*- coding:utf-8 -*-
import pandas as pd
import sys,os,glob,re
import argparse
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-gene",dest="gene")
parser.add_argument("-none",action="store_true")
args=parser.parse_args()

genelist = ["DDR2","MTOR","MYC","FGFR2","AKT1","FGFR1","FGFR3","FLT3","CDK6","CDKN2A","HER2","EGFR","ERBB2","CDK4","PDGFRA","APC","MET","ATM","KDR"]
genes=[]
for inputgene in args.gene.split(","):
	if inputgene in genelist:
		print u"存在gene：",inputgene
	else:
		genes.append(inputgene)
gene = genes
print gene
cnv_add=pd.DataFrame({
	"gene":gene,
	"id":".",
	"panal":"518gene_3M"

},columns=["gene","id","panal"])
# cnv_add={
# 	"gene":gene,
# 	"id":".",
# 	"panal":"518_3M"
#
# }
print cnv_add

# all_file=glob.glob("{}/*.txt.all".format(os.getcwd()))[0]

# cnv_all=pd.read_table(all_file)
def merge():
	cnv_list=pd.read_table("/gpfs/software/pipeline/pipeline/database/CNV/genels_tag/518gene_3M.tag.CNV.ls",header=None,names=["gene","id","panal"])

	cnv_list.columns=["gene","id","panal"]
	print cnv_list
	if args.none:
		df=cnv_list
	else:
		df=pd.concat([cnv_list,cnv_add])

	print df
	os.system("mkdir cnv_result")

	df.to_csv("{}/518gene_3MtoAgilent_60M.tag.CNV.ls".format(os.getcwd()),sep="\t",index=False,header=False)
	# # print ds.indicator_conlumn.value_counts()
	sh=pd.read_table("cnvkit_filter_job.sh",header=None)
	print sh
	import re
	sh1=re.sub("\S+Agilent60M_gene_list","{}/518gene_3MtoAgilent_60M.tag.CNV.ls".format(os.getcwd()),sh[0][4])
	sh2=re.sub("--outdir \S+","--outdir {}/cnv_result".format(os.getcwd()),sh1)
	print sh2

	cmd = sh2
	return cmd,df

def MergePicfile():
	global dfs
	xlsdir = glob.glob("./*_ffpe_ffpe.CNVkit.paired.txt.for_pic.xls")[0]
	name = xlsdir.split("\/")[-1]
	df = pd.read_excel(xlsdir)
	genelist = ["DDR2","MTOR","MYC","FGFR2","AKT1","FGFR1","FGFR3","FLT3","CDK6","CDKN2A","HER2","EGFR","ERBB2","CDK4","PDGFRA","APC","MET","ATM","KDR"]

	dk = df[(df[u"基因"] =="DDR2")|(df[u"基因"] =="KIT")|(df[u"基因"] =="MTOR")|(df[u"基因"] =="MYC")|(df[u"基因"] =="FGFR2")|(df[u"基因"] =="AKT1")|(df[u"基因"] =="FGFR1")|(df[u"基因"] =="FGFR3")|(df[u"基因"] =="FLT3")|(df[u"基因"] =="CDK6")|(df[u"基因"] =="CDKN2A")|(df[u"基因"] =="HER2")|(df[u"基因"] =="EGFR")|(df[u"基因"] =="ERBB2")|(df[u"基因"] =="CDK4")|(df[u"基因"] =="PDGFRA")|(df[u"基因"] =="APC")|(df[u"基因"] =="MET")|(df[u"基因"] =="ATM")|(df[u"基因"] =="KDR")]

	if args.gene:
		for i in gene:
			dk = pd.concat([dk,df[(df[u"基因"] ==i)]],axis=0)
			print df[(df[u"基因"] ==i)]
	# print dk[u"拷贝率是否变异"]
	pgenes=[]
	for pgene ,indexs in zip(dk[dk[u"拷贝率是否变异"]==u"是"][u"基因"],dk[dk[u"拷贝率是否变异"]==u"是"][u"基因"].index):
		if pgene in args.gene.split(","):
			# pgenes.append(pgene)
			pass
		else:
			pgenes.append(pgene)
			dk = dk.drop(indexs)
			# dfs = dfs.drop(indexs)
	print gene,pgenes
	print dk
	dfs=dfs[-dfs.gene.isin(pgenes)]
	dk.to_excel(glob.glob("./cnv_result/")[0]+name,index=False)

	# print genelist #in df[u"基因"]

cmd ,dfs= merge()

sampleid = re.findall("--sample_id\s\S+",cmd)[0].split(" ")[1]

workfile = glob.glob("/gpfs/share/bioinfo/商业样本/WES/*/*/{sampleid}*".format(**locals()))[0]
workfile = re.sub("QC_REPORT","CNV",workfile)
fatherpath = os.path.abspath(os.path.dirname(workfile)+os.path.sep+".")

cmd += "\n\ncp ./cnv_result/* {fatherpath} \n\nchmod -R 775 {fatherpath}".format(**locals())
# cmd+= "/gpfs/home/songtingrui/src_bio/cnv/nv/getData.sh  /gpfs/project/RD_bioinfor/WES_clinical/20180926_WES/wangtong_ffpe/cnv/tumor/cnvkit/it/PM2018091523.cnr /gr /gpfs/project/RD_bioinfor/WES_clinical/20180926_WES/wangtong_ffpe/cnv/normal/cnvkit/it/PM2018091507.cnr /gr /gpfs/project/RD_bioinfor/WES_clinical/20180926_WES/wangtong_ffpe/cnv/cnv_result cnvresul"
MergePicfile()
dfs.to_csv("{}/518gene_3MtoAgilent_60M.tag.CNV.ls".format(os.getcwd()),sep="\t",index=False,header=False)
files = open("cnvmergepic.sh","w")
files.write(cmd)

os.system(cmd)