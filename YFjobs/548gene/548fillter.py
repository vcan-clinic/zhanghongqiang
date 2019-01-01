import pandas as pd
import sys,re,os,glob
# vipdata = pd.read_table("/gpfs/home/zhaohongqiang/project/548gene/result/sunqixiang_ffpe_vip.xls",sep="\t")
annopath = glob.glob("/gpfs/project/RD_bioinfor/KY_project/RD_548gene/20180922_548gene/*/*/annotation/")

for i in annopath:
	print i
	name = re.findall("\w+_ffpe|\w+_cf", i)[0]
	cmd = "python /gpfs/home/zhaohongqiang/project/548gene/annovar_annotation_filter.py {i}*3tool.vcf.gz {i}*3tool.hg19_multianno.txt ~/scripts/mk_report/filter_config.yaml /gpfs/home/zhaohongqiang/project/548gene/result/{name}_vip.xls".format(**locals())
	try:
		os.system(cmd)
	except:
		print "{i} no work".format(**locals())
def toexceldata(vipdata):
	name = re.findall("\w+_ffpe|\w+_cf", vipdir)[0]
	vipresult = vipdata[
		((vipdata["normal_ratio"] == 0) & (vipdata["tumor_ratio"] >= 1)) |
		((vipdata["normal_ratio"] <= 1) & (vipdata["tumor_ratio"] / vipdata["normal_ratio"] >= 10))
		]
	# vipresult.to_csv("/gpfs/home/zhaohongqiang/project/548gene/result/test.csv")
	vipresult2 = vipresult[
		[u'chr', u'pos', u'ref', u'alt', u'Gene.refGene', u'exon', u'tumor_ratio', u'normal_ratio', u'transcripts',
		 u'CDS_Mutation', u'AA_Mutation', u'alt_depth', u'all_depth', u'Func.refGene', u'ExonicFunc.refGene',u"GeneDetail.refGene",
		 u'cosmic70']
	]
	vipresult2.drop_duplicates([u'chr',
								u'pos', u'ref', u'alt', ]).to_csv(
		"/gpfs/home/zhaohongqiang/project/548gene/result/{name}.548gene.csv".format(**locals()), index=False)

def toexceldata_cf(vipdata):
	name = re.findall("\w+_ffpe|\w+_cf", vipdir)[0]
	vipresult = vipdata[
		((vipdata["normal_ratio"] == 0) & (vipdata["tumor_ratio"] >= 0.5)) |
		((vipdata["normal_ratio"] <= 1) & (vipdata["tumor_ratio"] / vipdata["normal_ratio"] >= 10))
		]
	# vipresult.to_csv("/gpfs/home/zhaohongqiang/project/548gene/result/test.csv")
	vipresult2 = vipresult[
		[u'chr', u'pos', u'ref', u'alt', u'Gene.refGene', u'exon', u'tumor_ratio', u'normal_ratio', u'transcripts',
		 u'CDS_Mutation', u'AA_Mutation', u'alt_depth', u'all_depth', u'Func.refGene', u'ExonicFunc.refGene',u"GeneDetail.refGene",
		 u'cosmic70']
	]
	vipresult2.drop_duplicates([u'chr',
								u'pos', u'ref', u'alt', ]).to_csv(
		"/gpfs/home/zhaohongqiang/project/548gene/result/{name}.548gene.csv".format(**locals()), index=False)
vippath = glob.glob("/gpfs/home/zhaohongqiang/project/548gene/result/*vip.xls")
for vipdir in vippath:

	name = re.findall("\w+_ffpe|\w+_cf", vipdir)[0]
	if "_ffpe" in vipdir:
		vipdata=pd.read_table(vipdir,sep="\t")
		toexceldata(vipdata)
	if "_cf" in vipdir:
		vipdata = pd.read_table(vipdir, sep="\t")
		toexceldata_cf(vipdata)