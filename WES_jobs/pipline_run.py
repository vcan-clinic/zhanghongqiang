#-*- coding:utf-8 -*-
import os,re,sys,time
import paramiko
import pandas as pd
import argparse
from collections import defaultdict
import glob
import hashlib
def makehash():
    return defaultdict(makehash)
Dicter=makehash()
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-company",dest="company",help=u"请首先输入测序公司novo / fsz/no" ,default="novo")
parser.add_argument("-fp",action="store_true",help=u"``````````````````````````开启fastp ```````````````````` \n例：")
parser.add_argument("-ffpePM",dest="fpPM",help=u"ffpePM号")
parser.add_argument("-wcPM",dest="wcPM",help=u"WCPM号")
parser.add_argument("-fpinput",dest="fpinput",help="input")
# parser.add_argument("-company",dest="company",help=u"测序公司")
# parser.add_argument("-fpoutput",dest="fpoutput",help="output")



parser.add_argument("-run",action="store_true",help=u"````````````````开启run.pipeline``````````````````````\n例：",)
parser.add_argument("-o","--outfile",dest="outfile",help="outfile")
parser.add_argument("-s","--samplesheet",dest="samplesheet",default="SampleSheet.csv",help="samplesheet")
parser.add_argument('-fq',"--fastq",dest='FastqDir',help=u"fastaqdir")
parser.add_argument("-panal",dest="panal",default="Agilent_60M",help="input Agilent_60M,518gene_3M,548gene_3.2M")
parser.add_argument("-trrads",dest="trrads",default=8,help=u"核数")
parser.add_argument("-analysis",dest="analysis",default="qc,snv,cnv,fusion",help=u"分析项目")
parser.add_argument("-notpair",dest="notpair",action="store_true",help=u"是否为pair，默认是")

# parser.add_argument("-")`
parser.add_argument("-send",action="store_true",help=u"`````````````````````上传共享盘``````````````````````````\n例：")
parser.add_argument("-sample",dest="sample",help=u"文件夹")



parser.add_argument('-down',action="store_true",help=u"````````````````````开启下载，默认关闭`````````````````````\n例：nohup "  )
parser.add_argument("-i", '--ip', dest="host", default="36.110.13.43",
					  help="sftp ip ,default 36.110.13.43 ")
parser.add_argument("-p", '--port', dest="port", default=11211,
					  help="port , default 11211")
parser.add_argument("-u", '--user', dest="user", default="ZYXHSO37",
					  help="user , default ZYXHSO37")
parser.add_argument("-w", '--password', dest="password", default="ZYXH$$fsz1234",
					  help="password, default ZYXH$$fsz1234")
parser.add_argument("-d", '--localdir', dest="localdir", default=True,
					  help=u"the dir rawdata save in your machine本地地址")
parser.add_argument("-r", '--remotedir', dest="remotedir", default=True,
					  help=u"the dir in server machine ftp服务器地址")
# parser.add_option("-fq",'--fastq',dest='fastq')
args = parser.parse_args()
if args.company=="no":
	args.company ="fsz"
#############读取samplesheet#################
def read_samplesheet(sample_sheet=args.samplesheet):
	# #
	sampile_sheet_dict=makehash()
	# sampleSHEET=pd.read_csv(sample_sheet,header=20)
	sampledir=os.path.abspath(args.samplesheet)
	# PMid=sampleSHEET[19:]["Unnamed: 1"]
	# panal= sampleSHEET[19:]["Unnamed: 10"]
	# print sampleSHEET
	if args.company=="novo":
		sampleSHEET = pd.read_csv(sample_sheet, header=20)
		for PM,name,panal in zip(sampleSHEET["Sample_Name"],sampleSHEET["Description"],sampleSHEET["panel"]):
			# sampile_sheet_dict["PM"]=PM

			sampile_sheet_dict[name]["PM"]=PM
			sampile_sheet_dict[name]["panal"]=panal
		print sampile_sheet_dict
		return sampledir, sampile_sheet_dict
	else:
		sampleSHEET = pd.read_csv(sample_sheet, header=19)
		for PM,name,panal in zip(sampleSHEET["Sample_Name"],sampleSHEET["Description"],sampleSHEET["panel"]):
			# sampile_sheet_dict["PM"]=PM

			sampile_sheet_dict[name]["PM"]=PM
			sampile_sheet_dict[name]["panal"]=panal
		# print sampile_sheet_dict
		print sampile_sheet_dict
		return sampledir,sampile_sheet_dict





def run_pipline(fastqdir=args.FastqDir):
	ctrl={}
	fastq_dir=os.path.abspath(fastqdir)
	to_pipline=open("run.pipeline.sh","w")
	project_dir=os.getcwd()
	sampledir, sampledict = read_samplesheet()
	for name,sample in sampledict.items():
		print "name:",name , "	panal:",sample.values()[0] ,"	PM:",sample.values()[1] ,"fastq:",fastq_dir
	if args.notpair:
		ctrl["caller"]="vardict,vardictjava,varscan2,freebayes"
	else:
		ctrl["caller"]="mutect2,varscan2_pair,vardictjava_pair"
	if args.panal=="Agilent_60M" :
		analysis = args.analysis
		ctrl["pipline"] = "/gpfs/home/zhaohongqiang/pipeline/pipeline_v3.5_clinical_exome/run_pipeline.py"
		ctrl["config_file"]="/gpfs/home/zhaohongqiang/pipeline/pipeline_v3.5_clinical_exome/conf_lsf_local_yanfa.exon.yaml"
		ctrl["wai_bao"] = "\\\n--wai_bao 2000 \\\n"
	if args.panal in"518gene_3M":
		ctrl["pipline"] ="/gpfs/software/pipeline/pipeline/pipeline/run_pipeline.py"
		ctrl["config_file"]="/gpfs/software/pipeline/pipeline/pipeline/conf_lsf_local.yaml"
		ctrl["wai_bao"]="\\\n"
		analysis=args.analysis+"contamination"
	if args.panal=="548gene_3.2M" or args.panal in ['31gene_IG',"P13","P2","P1_2","IDT_pool1234","21gene_135k",'59gene_177K', '59gene_26K', '13gene_28K', '518gene_3M', '59gene_188K',
                                         '59gene_115K', '14gene_prostate', '15gene_liver', '18gene_gastric', '30gene_gyn',
                                         '38gene_breast','breast_148K', 'prostate_148K', 'liver_114K', 'gastric_114K', 'gyn_101K']:
		analysis = args.analysis
		ctrl["pipline"]="/gpfs/home/zhaohongqiang/pipeline/pipeline_v3.5_548gene/run_pipeline.py"
		ctrl["config_file"]="/gpfs/home/zhaohongqiang/pipeline/pipeline_v3.5_548gene/conf_lsf_local_yanfa.yaml"
		ctrl["wai_bao"]="\\\n --wai_bao 2000  \\\n"
	# print ctrl
	print "caller:",ctrl["caller"],"\npipline:",ctrl["pipline"],"\nproject_dir:",project_dir
	run_pipline_write="""#!/bin/bash\n/gpfs/software/bcbio/bcbio/anaconda/bin/python {pipline} \\\n--project_dir {project_dir} \\\n--fastq_dir {fq_dir} \\\n--sample_sheet {samplesheet_dir} \\\n--config_file {config_file} \\\n--num_threads {trrads} \\\n--core_node 1 \\\n--scheduler lsf \\\n--analysis {analysis} \\\n--caller {caller} \\\n--chip {chip} {wai_bao}--to_db no""".format(pipline=ctrl["pipline"],project_dir=project_dir,fq_dir=fastq_dir,samplesheet_dir=sampledir,config_file=ctrl["config_file"],trrads=args.trrads,analysis=analysis,caller=ctrl["caller"],chip=args.panal,wai_bao=ctrl["wai_bao"])
	to_pipline.write(run_pipline_write)
	to_pipline.close()


def Fastp_job():
	args.samplesheet = os.path.abspath(os.path.join(args.fpinput, ".."))
	sampledir, sampledict = read_samplesheet()
	global PM, wcPM
	for name, sample in sampledict.items():
		if "FFPE" in name or "ffpe" in name:
			PM = sample["PM"]
		if "WC" in name or "wc" in name:
			wcPM = sample["PM"]
	fastq_input=os.path.abspath(args.fpinput)
	fastq_output=os.path.abspath(args.fpinput)
	os.system("mkdir {fastq_input}/fastp_splitfastq/".format(**locals()))
	def findfastq(company):
		if company=="novo":
			fastq_R1_dir=glob.glob("{fastq_input}/{PM}*L*_1*".format(fastq_input=fastq_input,PM=PM))[0]
			# print fastq_R1_dir
			fastq_R2_dir=glob.glob("{fastq_input}/{PM}*L*_2*".format(fastq_input=fastq_input,PM=PM))[0]
			return fastq_R1_dir, fastq_R2_dir
		elif company=="fsz":
			fastq_R1_dir = glob.glob("{fastq_input}/{PM}*R1*".format(fastq_input=fastq_input, PM=PM))[0]
			# print fastq_R1_dir
			fastq_R2_dir = glob.glob("{fastq_input}/{PM}*R2*".format(fastq_input=fastq_input, PM=PM))[0]
			return 	fastq_R1_dir,fastq_R2_dir
	fastq_R1_dir,fastq_R2_dir = findfastq(company=args.company)
	fastq_dir_R1_out=args.fpinput+"/fastp_splitfastq/{dir}".format(PM=PM,dir=fastq_R1_dir.split("/")[-1])
	fastq_dir_R2_out = args.fpinput + "/fastp_splitfastq/{dir}".format(PM=PM, dir=fastq_R2_dir.split("/")[-1])
	fastq_dir_out=args.fpinput+"/fastp_splitfastq"
	html_dir=args.fpinput+"/fastp_splitfastq/splitresul{PM}.html".format(PM=PM,dir=fastq_R1_dir.split("/")[-1])
	jeson_dir=args.fpinput+"/fastp_splitfastq/splitresul{PM}.jeson".format(PM=PM,dir=fastq_R1_dir.split("/")[-1])

	fastp_job = open("fastp{PM}_jobs.sh".format(**locals()),"w")
	cmd = "/gpfs/home/liufx/pipeline/tools/fastp -i {fastq_R1_dir} -I {fastq_R2_dir} -s 6 -o {fastq_dir_R1_out} -O {fastq_dir_R2_out} -A -Q -G -h {html_dir} -j {jeson_dir} >run.log".format(**locals())
	cmd +="\ncp -d {wcPM}* {fastq_dir_out}".format(**locals())
	fastp_job.write(cmd)
	# os.system(cmd)
	fastp_job.close()

	# os.system("cp -d {wcPM}* {fastq_dir_out}".format(**locals()))
	for i in  glob.glob("{fastq_dir_out}/000*{PM}*".format(**locals())):

		# if  re.findall(r"000\S+",i)[0]:

		if args.company=="novo":
			old_name = re.findall(r"PM\S+", re.findall(r"000\S+", i)[0])[0].replace("_L\d", "_L{}".format(
				re.findall(r"000\S+", i)[0].split(".")[0]))
			# print  re.search("(\d+).(RD\d+|PM\d+).*_L(\d+)_(\d+).fq.gz", i).groups()
			nm,rd, l, r = re.search("(\d+).(RD\d+|PM\d+).*_L(\d+)_(\d+).fq.gz", i).groups()
			new_name = "{0}_S4_L{1}_{2}.fastq.gz".format(rd, l+nm, r)
			new_name_dir=os.path.join(fastq_dir_out,new_name)
			try:
				os.symlink(i,new_name_dir)
			except:
				print "File exists"
		else:
			os.system ("mv {fastq_dir_out}/{bl} {new}/{out}".format(fastq_dir_out=fastq_dir_out,new=fastq_dir_out,bl=re.findall(r"000\S+",i)[0]  , out= re.findall(r"PM\S+",re.findall(r"000\S+",i)[0])[0].replace("L00","L{}".format(re.findall(r"000\S+",i)[0].split(".")[0])) ))
	if args.company == "novo":
		for i in glob.glob("{fastq_dir_out}/{wcPM}*".format(**locals())):
			rd, l, r = re.search("(RD\d+|PM\d+).*_L(\d+)_(\d+).fq.gz", os.path.split(i)[1]).groups()
			new_name = "{0}_S4_L{1:0>3}_{2}.fastq.gz".format(rd, l, r)
			new_name_dir = os.path.join(fastq_dir_out, new_name)
			print i,new_name_dir
			try:
				os.symlink(i,new_name_dir)
			except:
				print "File exists"
def make_518resul():
	pass













########################FTP文件下载##############################
def DOWN_sftp(host, port, username, password, local, remotes):
		sf = paramiko.Transport((host, int(port)))
		sf.connect(username=username, password=password)
		sftp = paramiko.SFTPClient.from_transport(sf)
		local_dir_files = os.listdir(local)
		num = 1
		try:
			for remote in remotes:
				if os.path.isdir(local):  # 判断本地参数是目录还是文件
					for f in sftp.listdir(remote):  # 遍历远程目录
						# print os.path.join(remote,f)
						if f == "nohup.out": continue
						if f.endswith("txt"):
							sftp.get(os.path.join(remote, f), os.path.join(local, "md5_%s.txt" % num))
							num += 1
							continue
						if f in local_dir_files: continue
						sftp.get(os.path.join(remote, f), os.path.join(local, f))  # 下载目录中文件
					# print(f)
				else:
					# sftp.get(remote,local)#下载文件
					raise Exception("local is not dir")
		except Exception, e:
			print('download exception:', e)
		sf.close()


def md5sum(fname):
    """ 计算文件的MD5值
    """
    def read_chunks(fh):
        fh.seek(0)
        chunk = fh.read(8096)
        while chunk:
            yield chunk
            chunk = fh.read(8096)
        else: #最后要将游标放回文件开头
            fh.seek(0)
    m = hashlib.md5()
    if isinstance(fname, basestring) \
            and os.path.exists(fname):
        with open(fname, "rb") as fh:
            for chunk in read_chunks(fh):
                m.update(chunk)
    #上传的文件缓存 或 已打开的文件流
    elif fname.__class__.__name__ in ["StringIO", "StringO"] \
            or isinstance(fname, file):
        for chunk in read_chunks(fname):
            m.update(chunk)
    else:
        return ""
    return m.hexdigest()


def down_data():
	args = parser.parse_args()
	host = args.host
	port = args.port
	USER = args.user
	PASSWORD = args.password
	local = args.localdir
	remot1 = args.remotedir.strip().split(",")

	# print host,port,USER,PASSWORD,local,remot1

	f=open("local_md5.txt","w")

	DOWN_sftp(host, port, USER, PASSWORD, local, remot1)
	for i in glob.glob(local+"/*.fastq.gz"):
		f.write(md5sum(i)+"\t"+i)
	f.close()
	os.system("touch %s" % os.path.join(local, "RTARead3Complete.txt"))

def SendFile():

	gxfile_wes = "/gpfs/share/bioinfo/商业样本/WES/"
	name = args.sample
	samplie = glob.glob(os.path.abspath(args.sample))[0]  # 数据文件的绝对路径
	time_file = time.strftime('%Y%m%d', time.localtime(time.time()))
	htmlfile = glob.glob(samplie+"/qc_report/*_ffpe.qc.html")[0]
	Allfile = glob.glob(samplie+"/annotation/*_somatic.xls")[0]
	samplecold = re.findall(r"WES_\w+(?=_somatic)",Allfile.split("annotation")[-1])[0]
	PM = re.findall(r"PM\w+(_ffpe)",htmlfile)[0]
	directory = "{samplie}/{time_file}directory".format(**locals())
	os.system("mkdir {directory}".format(**locals()))
	def ColuMns():
		somatic_col = [u'Samplecode', u'taged', u'Class', u'Reported', u'Score', u'532GENE',
				   u'Software', u'GENE', u'CDS_Mutation', u'AA_Mutation', u'Exon',
				   u'tumor_ratio%', u'normal_ratio%', u'Leve', u'QC', u'Likely_germline',
				   u'NM', u'status', u'target_chemo_mark', u'COSMIC',
				   u'COSMIC_somatic_status', u'ExonicFunc', u'alt_depth', u'all_depth',
				   u'control_alt_depth', u'control_all_depth',
				   u'tumor_alt_forward_reverse', u'normal_alt_forward_reverse',
				   u'homopolymer', u'dbsnp', u'1000G_ALL', u'1000G_EAS', u'ExAC_Freq',
				   u'ExAC_EAS', u'ESP6500si_ALL', u'ClinVar_SIG', u'ClinVar_DIS',
				   u'Pathogenicity', u'soft_count', u'Freebayes_num', u'Vardict_num',
				   u'Mutation_db', u'DISEASE', u'FDACFDA1', u'FDACFDA2', u'DRUGS',
				   u'INFORMATION', u'COSMIC_ID', u'chr', u'pos', u'gene_symbol', u'ref',
				   u'alt', u'SIFT_pred', u'MutationAssessor_pred', u'Polyphen2_HDIV_pred',
				   u'Polyphen2_HVAR_pred', u'LRT_pred', u'MutationTaster_pred',
				   u'FATHMM_pred', u'PROVEAN_pred', u'fathmm-MKL_coding_pred',
				   u'MetaSVM_pred', u'MetaLR_pred', u'是否报出', u'是否有药', u'功能', u'批准疾病',
				   u'用于批准疾病的敏感用药', u'用于其他疾病的敏感用药', u'临床试验用药', u'其他耐药信息', u'级别', u'注释',
				   u'Chr', u'Position', u'Ref', u'A', u'T', u'C', u'G', u'N']
		germline_col = [u'Samplecode', u'taged', u'Class', u'Reported', u'Score', u'532GENE',
					 u'Software', u'GENE', u'CDS_Mutation', u'AA_Mutation', u'Exon',
					 u'normal_ratio%', u'tumor_ratio%', u'Leve', u'QC', u'Likely_germline',
					 u'NM', u'status', u'target_chemo_mark', u'COSMIC',
					 u'COSMIC_somatic_status', u'ExonicFunc', u'alt_depth', u'all_depth',
					 u'control_alt_depth', u'control_all_depth',
					 u'tumor_alt_forward_reverse', u'normal_alt_forward_reverse',
					 u'homopolymer', u'dbsnp', u'1000G_ALL', u'1000G_EAS', u'ExAC_Freq',
					 u'ExAC_EAS', u'ESP6500si_ALL', u'ClinVar_SIG', u'ClinVar_DIS',
					 u'Pathogenicity', u'soft_count', u'Freebayes_num', u'Vardict_num',
					 u'Mutation_db', u'DISEASE', u'FDACFDA1', u'FDACFDA2', u'DRUGS',
					 u'INFORMATION', u'COSMIC_ID', u'chr', u'pos', u'gene_symbol', u'ref',
					 u'alt', u'SIFT_pred', u'MutationAssessor_pred', u'Polyphen2_HDIV_pred',
					 u'Polyphen2_HVAR_pred', u'LRT_pred', u'MutationTaster_pred',
					 u'FATHMM_pred', u'PROVEAN_pred', u'fathmm-MKL_coding_pred',
					 u'MetaSVM_pred', u'MetaLR_pred', u'临床意义', u'注释']
		cnv_col = [u"样本编号",u"基因",u"染色体",u"拷贝率",u"拷贝率是否变异"]
		return somatic_col,germline_col,cnv_col

	somatic_col, germline_col,cnv_col=ColuMns()
####wason文件####
	cmd = "mkdir {gxfile_wes}/WATSON/{time_file}/ ".format(**locals())
	cmd += "\n\n  mkdir {gxfile_wes}/WATSON/{time_file}/{name}".format(**locals())
	anno = glob.glob(samplie + "/annotation/*.sofstric.watson.vcf.hg19_multianno.txt")[0]
	_518dir = re.sub("_WES","_518",samplie)
	wason = glob.glob(samplie+"/*.collection/")[0] #wason文件地址
	cmd += "\ncp -r {anno} {gxfile_wes}/WATSON/{time_file}/{name}".format(**locals())
	cmd += "\ncp -r {wason}/* {gxfile_wes}/WATSON/{time_file}/{name}".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/WATSON/{time_file}".format(**locals())




####化疗####
	os.system("mkdir {directory}/chemo".format(**locals()))
	chemo =	glob.glob(samplie+"/annotation/*_chemo.xls")[0] #化疗文件地址
	chemofile = pd.read_excel(chemo)
	# print chemofile
	os.system("cp {chemo} {directory}/chemo".format(**locals()))
	chemodel20 = glob.glob("{directory}/chemo/*".format(**locals()))[0]
	chemofile2 = chemofile[(chemofile["Depth"]>=20)&(chemofile["Depth"]== chemofile["Depth"])]
	print chemofile2
	chemofile2.to_excel(chemodel20,index=False)
	# print chemofile2
	cmd += "\n\nmkdir {gxfile_wes}/CHEMO/{time_file}\ncp {chemodel20} {gxfile_wes}CHEMO/{time_file}".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/CHEMO/{time_file}".format(**locals())





####mutation#####
	make_weidian =os.path.abspath(os.path.dirname(samplie))
	# print make_weidian
	result ="/gpfs/software/anaconda2.4.3/bin/python   /gpfs/software/pipeline/pipeline/clinical_examination_script/Stat_result/make_result2.py -l -e -t --project_dir {make_weidian} --out_dir {make_weidian} ".format(**locals())
	if glob.glob(make_weidian+"/20*"):
		pass
	else:
		os.system(result)
	# weidian = glob.glob(make_weidian+"/20*/*/weidian")[0]
	# # result1 = " python /gpfs/home/shizhenyu/scripts/snv_filter.py {weidian} ".format(**locals())
	# # os.system(result1)
	# # mutation = glob.glob(samplie+"/annotation/")[0]
	# anno = glob.glob(samplie+"/annotation/*.sofstric.watson.vcf.hg19_multianno.txt")[0]
	# cmd += "\nmkdir {gxfile_wes}/MUTATION/{time_file} \n cp -r {weidian}/{name} {gxfile_wes}/MUTATION/{time_file}/".format(**locals())
	# cmd += "\ncp -r {anno} {gxfile_wes}/MUTATION/{time_file}/".format(**locals())
	# cmd += "\n chmod -R 775 {gxfile_wes}/MUTATION/{time_file}".format(**locals())





####cnv#####
	cnv = glob.glob(samplie + "/cnv/*paired.txt.for_pic.xls")[0]
	cnv_html = glob.glob(samplie + "/cnv/*_ffpe.html")[0]
	cnv_all = glob.glob(samplie + "/cnv/*_ffpe_ffpe.CNVkit.paired.txt.all")[0]
	# print cnv
	cmd +=	"\n\nmkdir {gxfile_wes}/CNV/{time_file}/ ".format(**locals())
	cmd +=	"\n\ncp -rd {cnv} {gxfile_wes}/CNV/{time_file}/".format(
		**locals())
	cmd +=	"\n\ncp -rd {cnv_html} {gxfile_wes}/CNV/{time_file}/".format(**locals())
	cmd += "\n\ncp -rd {cnv_all} {gxfile_wes}/CNV/{time_file}/".format(**locals())
	cmd +=	"\nchmod -R 775 {gxfile_wes}/CNV/{time_file}".format(**locals())
	# print [make_weidian],weidian





####QC#####
	QC = wason+"/QC"
	cmd += "\n\nmkdir {gxfile_wes}/QC/{time_file} && mkdir {gxfile_wes}/QC/{time_file}/{name} \n cp -rd {QC} {gxfile_wes}/QC/{time_file}/{name}".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/QC/{time_file}".format(**locals())





####QC_report###
	collection_log_data =glob.glob(samplie+"/collection.log.txt")[0]
	# pd.read_csv(collection_log_data,header=None,sep="\n")[0][1].split(",")
	os.system("mkdir {directory}/qc_report".format(**locals()))
	QC_Report_file = glob.glob("{directory}/qc_report".format(**locals()))[0]
	# print QC_Report_file
	collection_log_list = pd.read_csv(collection_log_data, header=None, sep="\n")[0][1].split(",")
	qc_data = pd.read_csv(collection_log_data, header=None, sep="\n")[0][1]
	# print qc_data
	if float(re.findall(r"""Raw Q30 rate\S+\s\w+\W\w+""", qc_data)[0].split(" ")[-1]) > 80:
		RQ = u"合格"
	else:
		RQ = u"不合格"
	if float(re.findall(r"""Coverage\s\S+\s\w+\W\w+""", qc_data)[0].split(" ")[-1]) > 85:
		fm = u"合格"
	else:
		fm = u"不合格"
	if RQ == fm == u"合格":
		ql = u"合格"
	else:
		q1 = u"不合格"
	qcdata = {
		"Samplecode":samplecold,
		# "QC": collection_log_list[0].split(" ")[1],
		"Average depth(rmdup)": re.findall(r"""Average depth\S+\s\w+\W\w+""", qc_data)[0].split(" ")[-1],
		"Coverage": fm,
		"Raw Q30 rate(%)": RQ,
		"Fraction of mapped bases(%)": re.findall(r"""Fraction of mapped bases\S+\s\w+\W\w+""", qc_data)[0].split(" ")[-1],
		"Quality": ql
	}
	# print qcdata
	QC_Report = pd.DataFrame(qcdata,index=[0],columns=["Samplecode","Raw Q30 rate(%)","Fraction of mapped bases(%)","Average depth(rmdup)","Coverage","Quality"])
	QC_Report.to_excel(QC_Report_file+"/"+samplecold+".xls",index=None)
	qc_report= glob.glob(QC_Report_file+"/"+samplecold+".xls")[0]
	cmd += "\nmkdir {gxfile_wes}/QC_REPORT/{time_file}/ \n cp -rd {qc_report} {gxfile_wes}/QC_REPORT/{time_file}/".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/QC_REPORT/{time_file}".format(**locals())




####TMB###
	cmd +="\nmkdir {gxfile_wes}/TMB/{time_file}/ ".format(**locals())
	# print PM
	tmb= glob.glob(samplie+"/annotation/*_TMB.xls".format(**locals()))[0]
	# print tmb
	cmd += "\ncp -r {tmb} {gxfile_wes}/TMB/{time_file}/".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/TMB/{time_file}".format(**locals())




####WES_mutation####
	read_genetic = glob.glob(samplie+"/annotation/READ*genetic.xls")[0]
	read_somatic = glob.glob(samplie+"/annotation/READ*somatic.xls")[0]
	ommifile = glob.glob(samplie+"/annotation/OMMI_*_genetic.xls")[0]
	cmd += "\nmkdir {gxfile_wes}/WES_MUTATION/{time_file}/".format(**locals())
	cmd += "\nmkdir {gxfile_wes}/WES_MUTATION/{time_file}/{name}".format(**locals())
	cmd += "\ncp -r {read_genetic} {gxfile_wes}/WES_MUTATION/{time_file}/{name}".format(**locals())
	cmd += "\ncp -r {read_somatic} {gxfile_wes}/WES_MUTATION/{time_file}/{name}".format(**locals())
	cmd += "\ncp -r {ommifile} {gxfile_wes}/WES_MUTATION/{time_file}/{name}".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/WES_MUTATION/{time_file}".format(**locals())




#####all_mutation###
	all_genetic = glob.glob(samplie + "/annotation/ALL*genetic.xls")[0]
	all_somatic = glob.glob(samplie + "/annotation/ALL*somatic.xls")[0]
	cmd += "\nmkdir {gxfile_wes}/ALL_MUTATION/{time_file}/".format(**locals())
	cmd += "\nmkdir {gxfile_wes}/ALL_MUTATION/{time_file}/{name}".format(**locals())
	cmd += "\ncp -r {all_genetic} {gxfile_wes}/ALL_MUTATION/{time_file}/{name}".format(**locals())
	cmd += "\ncp -r {all_somatic} {gxfile_wes}/ALL_MUTATION/{time_file}/{name}".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/ALL_MUTATION/{time_file}".format(**locals())




#####DMMR########
	all_dmmr_result = glob.glob(samplie + "/annotation/ALL*_dmmr_result.xls")[0]
	all_dmmr_sample = glob.glob(samplie + "/annotation/ALL*_dmmr_sample.xls")[0]
	cmd += "\nmkdir {gxfile_wes}/DMMR/{time_file}/".format(**locals())
	cmd += "\nmkdir {gxfile_wes}/DMMR/{time_file}/{name}".format(**locals())
	cmd += "\ncp -r {all_dmmr_result} {gxfile_wes}/DMMR/{time_file}/{name}".format(**locals())
	cmd += "\ncp -r {all_dmmr_sample} {gxfile_wes}/DMMR/{time_file}/{name}".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/DMMR/{time_file}".format(**locals())


#####MSI########
	cmd += "\nmkdir {gxfile_wes}/MSI/{time_file}/ ".format(**locals())
	msi = glob.glob(samplie + "/annotation/*_MSI.xls".format(**locals()))[0]
	cmd += "\ncp -r {msi} {gxfile_wes}/MSI/{time_file}/".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/MSI/{time_file}".format(**locals())

#####报点文件####
	cmd += "\nmkdir {gxfile_wes}/需要报出的点/{time_file}".format(**locals())

	try:
		os.mkdir("{directory}/post".format(**locals()))
	except:
		print "文件夹已经存在"
	writer = pd.ExcelWriter('{directory}/post/{name}.xlsx'.format(**locals()))
	somdf = pd.DataFrame(columns=somatic_col)
	gerdf = pd.DataFrame(columns=germline_col)
	cnvdf = pd.DataFrame(columns=cnv_col)
	somdf.to_excel(writer,sheet_name="somatic",index=None,encoding='utf-8')
	gerdf.to_excel(writer,sheet_name="germline",index=None,encoding='utf-8')
	cnvdf.to_excel(writer,sheet_name="CNV",index=None,encoding='utf-8')
	cmd += "\ncp {directory}/post/* {gxfile_wes}/需要报出的点/{time_file}".format(**locals())
	cmd += "\nchmod -R 775 {gxfile_wes}/需要报出的点/{time_file}".format(**locals())
	writer.save()

	f=open("sendfile{time_file}.sh".format(**locals()),"w")
	f.write(cmd)
	# os.system("sh sendfile.sh ")





def main():
	run_pipline()




if __name__ == '__main__':

	if args.down:
		down_data()

	if args.fp:
		Fastp_job()
	if args.run:
		main()
	if args.send:
		SendFile()
	# sampledir, sampledict = read_samplesheet()
	# global PM,wcPM
	# for name, sample in sampledict.items():
	# 	if "FFPE"  in name or "ffpe" in name :
	#
	# 		PM = sample["PM"]
	# 	if "WC" in name or "wc" in name :
	# 		wcPM = sample["PM"]
	# print  os.path.abspath(os.path.join(args.fpinput, ".."))