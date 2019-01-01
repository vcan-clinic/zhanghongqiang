#-*- coding:utf-8 -*-
import os,sys,re,glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input1",dest="input1",help=u"输入snv基因")
parser.add_argument("-input2",dest="input2",help=u"输入cnv基因")
parser.add_argument("-name",dest="name",help=u"样本名字")
parser.add_argument("-sample",dest="sample",help=u"")
parser.add_argument("-path",dest="path",help=u"workdir")

args=parser.parse_args()
workdir = args.path
postdir = "/gpfs/share/bioinfo/商业样本/WES/需要报出的点/"
postpath = glob.glob(postdir+"/*/"+args.name+"*")[0]

if postpath:#判断cnv是否存在
	snv = []
	cnv = []
	if pd.read_excel(postpath,sheet_name="somatic")["GENE"].size>0 :
		for genesnvsomatic in pd.read_excel(postpath, sheet_name="somatic")["GENE"]:
			snv.append(genesnvsomatic)
	if pd.read_excel(postpath,sheet_name="germline")["GENE"].size>0:
		for genesnvgearmline in pd.read_excel(postpath, sheet_name="germline")["GENE"]:
			snv.append(genesnvgearmline)
	if pd.read_excel(postpath,sheet_name="CNV")[u"基因"].size>0 :
		print pd.read_excel(postpath,sheet_name="CNV")[u"基因"].size>0
		for genesnvgearmline in pd.read_excel(postpath, sheet_name="CNV")[u"基因"]:
			cnv.append(genesnvgearmline)
else:
	snv =args.input1
	cnv= args,input2

print "这是SNV：",snv
print "这是CNV:",cnv
mutfiles = glob.glob(workdir+"/annotation/*mut")[0]

cnv_all = glob.glob(workdir+"/cnv/*.all")[0]



def getgene(genelist,cril):
	# print  genelist
	# genedata = genelist.split(",")
	genedf = pd.DataFrame(genelist, columns=["gene"])["gene"]
	genedf.to_csv("{}/{}_mmapgenelist.txt".format(os.getcwd(), args.sample), sep="\t", index=False, header=False)
	if cril=="snv":
		genedf.to_csv("{}/{}_mmapsnvgenelist.txt".format(os.getcwd(), args.sample), sep="\t", index=False, header=False)
		snvgenefile = glob.glob(workdir + "/annotation/*mmapsnvgenelist.txt".format(cril))[0]
		return snvgenefile
	else:
		genedf.to_csv("{}/{}_mmapcnvgenelist.txt".format(os.getcwd(), args.sample), sep="\t", index=False, header=False)
		cnvgenefile = glob.glob(workdir + "/annotation/*mmapcnvgenelist.txt".format(cril))[0]
		return cnvgenefile
if snv:
	snvgenefile=getgene(snv,cril="snv")
	
else:
	snvgenefile = open("{}/{}_mmapsnvgenelist.txt".format(os.getcwd(), args.sample),"w")
	snvgenefile.close()
	snvgenefile = "{}/{}_mmapsnvgenelist.txt".format(os.getcwd(), args.sample)
if cnv:
	cnvgenefile=getgene(cnv,cril="cnv")
else:
	cnvgenefile = open("{}/{}_mmapcnvgenelist.txt".format(os.getcwd(), args.sample),"w")
	cnvgenefile.close()
	cnvgenefile = "{}/{}_mmapcnvgenelist.txt".format(os.getcwd(), args.sample)
sampleid = args.sample

runsh = open("{workdir}/annotation/circos_run_jobs.sh".format(**locals()),"w")
runjob = "/gpfs/home/zhaohongqiang/pipeline/pipe_bin/plot.cnv_snv.circos_v3 \
        -m {mutfiles} \
        -c {cnv_all} \
        -p {sampleid}_mmap \
        -g {snvgenefile} {cnvgenefile}".format(**locals())
runjob += "\n\nmv {workdir}/annotation/{sampleid}_mmap.circos.png {workdir}/annotation/{sampleid}_mmap.png".format(**locals())
# runsh.write(runjob)
# os.system(runjob)
GXfile = glob.glob("/gpfs/share/bioinfo/商业样本/WES/*/*/{sampleid}*".format(**locals()))[0]
workfile = re.sub("QC_REPORT|CNV","Mmap",GXfile)
print workfile
fatherpath = os.path.abspath(os.path.dirname(workfile)+os.path.sep+".")
print fatherpath
runjob += "\n\nmkdir {fatherpath}\n\ncp {workdir}/annotation/{sampleid}_mmap.png {fatherpath} \n\nchmod -R 775 {fatherpath}".format(**locals())
os.system(runjob)
runsh.write(runjob)
runsh.close()