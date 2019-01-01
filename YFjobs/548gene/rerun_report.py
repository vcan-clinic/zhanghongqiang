#coding=utf-8
import sys,re,glob,os
import collections as cl

def read_f(f):
	with open(f) as fo:
		return [x.strip().split("\t") for x in fo]
def read_run_list(run_list): #上机日期\t报告日期
	d = cl.defaultdict(dict)
	for line in read_f(run_list):
		RD,sample_code,detect_date,report_date = line
		d[detect_date][RD] = [sample_code,report_date]
	return d 
def find_aim(pd,out_path,info,scripts_path):
	fusion = glob.glob(os.path.join(pd,"fusion/tumor/lumpy","*tumor_fusion_decision.txt"))[0]
	RD = re.search("_(RD\S+)_tumor_fusion_decision.txt",fusion.split("/")[-1]).group(1)
	vcf = glob.glob(os.path.join(pd,"annotation","*reformed.*tool.vcf.gz"))[0]
	txt = glob.glob(os.path.join(pd,"annotation","reformed.*tool.hg19_multianno.txt"))[0]
	out_snv_report = os.path.join(out_path,RD+".report.xls")
	out_pdf_report = os.path.join(out_path,RD+".report.pdf")
	out_json_report = os.path.join(out_path,RD+".report.json")
	cmd_lines = []	
	print RD,info
	if RD not in info:
		print "#cant find {} in run_list".format(RD)
	else:
		sample_code,report_date = info[RD]
		cmd_lines.append("python {scripts_path}/annovar_annotation_filter.py {vcf} {txt} {scripts_path}/filter_config.yaml {out_snv_report}".format(**locals()))
		cmd_lines.append("python {scripts_path}/make_report.py -p {pd} -s {RD} -f {fusion} -r {out_snv_report} -c {scripts_path}/report.yaml -d {out_path}".format(**locals()))
		cmd_lines.append("python {scripts_path}/make_pdf_report.py {out_json_report} {out_pdf_report} {sample_code} {RD} {report_date}".format(**locals()))
	return cmd_lines
def find_sample(sample_list,run_list):
	run_list = read_run_list(run_list)
	cmd_lines = []
	for i in read_f(sample_list):
		sample = i[0]
		m = re.search("^(.*)/project_\d+samples_(\d+).*/13gene_2\dK/([^_]+)",sample)
		if m:
			d,date,s_name = m.groups()
			if date in run_list:
				out_dir = date+"_pdf_report"
				if not os.path.isdir(out_dir):
					os.mkdir(out_dir)	
				cmd_lines += find_aim(sample,os.path.abspath(out_dir),run_list[date],scripts_path)
		else:
			print sample
	with open("pdf_report.sh","w") as sh:
		print >>sh,"\n".join(cmd_lines)
scripts_path = "/gpfs/home/shizhenyu/scripts/mk_report"
find_sample(sys.argv[1],sys.argv[2]) #sample_list ,all_info.xls
