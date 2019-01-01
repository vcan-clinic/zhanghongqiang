#coding=utf-8
import sys
from qc_report import QC_report
# import annovar_annotation_filter
import yaml
import os
import glob
import re
import collections as cl
import argparse
import json 
def make_qc_report(project_dir,qc_cutoff,vip_sites,sample_code,results_dir):
	q = QC_report(project_dir,"",vip_sites)
	q.get_q30("tumor")
	qc_d,site_depth_d = q.get_maping_data()
	qc_judge=[]
	qc_rules = {"Median_depth": lambda x:float(qc_d[x]) >qc_cutoff["Median_depth"],
		    "hot_spot_min_depth" : lambda x:float(qc_d[x])>qc_cutoff["hot_spot_min_depth"],
		    "Raw_Q30 %" : lambda x:float(qc_d[x]) >qc_cutoff["Raw_Q30 %"],
		    "Fraction of mapped bases(%)" :lambda x:float(qc_d[x]) >qc_cutoff["Fraction of mapped bases(%)"],
		    "Total_bases" :lambda x:float(qc_d[x]) >qc_cutoff["Total_bases"]
			}
	report_out = open(os.path.join(results_dir,sample_code + "_qc_report.xls"),"w")
	header = sorted(qc_d.keys())
	print >>report_out,"qc_data\tvalue\tpass_or_fail"
	for i in header:
		target_judge = qc_rules.get(i,lambda x:"-")(i)
		if not target_judge:
			qc_judge.append(i)
			print >>report_out,"{0}\t{1}\t{2}".format(i,qc_d[i],"fail")
		else:
			print >>report_out,"{0}\t{1}\t{2}".format(i,qc_d[i],"pass")
	report_out.close()
	return qc_judge 


def make_snv_report(snv_txt,snv_sites,sample_code ,results_dir):
	snv_out = open(os.path.join(results_dir,sample_code + "_snv_report.xls"),"w")
	need_data=["tumor_ratio","alt_depth","all_depth"]
	print >>snv_out,"\t".join(["gene","cosmic_id","cds_mutation","aa_mutation","detected"]+need_data)
	if not snv_txt:
		print >>snv_out,"can not read snv report.xls."
		snv_out.close()
		return "can not find report.xls",0
	snvs = {}
	with open(snv_txt) as snv:
		header = snv.readline().strip().split("\t")
		for line in snv:
			d = dict(zip(header,line.strip().split("\t")))
			d["AA_Mutation"]=d["AA_Mutation"].split("(")[0]
			if d["tag"] == "unknown_19del" and ("19del" not in snvs):
				snvs["19del"] = d
			elif d["tag"] == "vip_19del": #vip 
				d["standard_info"] = ["19del({0})".format(x) for x in snv_sites["EGFR"]["19del"][(d["chr"],d["pos"],d["ref"],d["alt"])]]
				snvs["19del"] = d
			else:
				snvs[(d["chr"],d["pos"],d["ref"],d["alt"])] = d
	genes = sorted(snv_sites.keys())
#	sites = sorted(snv_sites.keys(),key=lambda x:(snv_sites[x][0],int(x[1]))) 
	snv_out_list=[]
	snv_detacted = 0
	for gene in genes:
		sites = sorted(snv_sites[gene].keys()) 
		for site in sites:
			if site in snvs:
				if site == "19del":
	                                site_info =[snvs[site][x] for x in  need_data]
					pipline_site_info = ["19del({0})".format(snvs[site][x].split(",")[0]) for x in ["cosmic70","CDS_Mutation","AA_Mutation"]]
					out_info = snvs[site].get("standard_info",pipline_site_info)
					print >>snv_out,"\t".join([gene]+out_info+["yes"]+site_info)
					snv_out_list.append([gene,"-","Exon 19 Deletion","Exon 19 Deletion"]+["yes"])
					snv_detacted += 1 
				else:
					site_info =[snvs[site][x] for x in  need_data]
					print >>snv_out,"\t".join([gene]+snv_sites[gene][site]+["yes"]+site_info)
					snv_out_list.append([gene]+snv_sites[gene][site]+["yes"])
					snv_detacted += 1
			else:
                                if site == "19del":
                                    print >>snv_out,"\t".join(["EGFR","-","Exon 19 Deletion","Exon 19 Deletion"]+["no"]+["-"]*len(need_data))
                                    snv_out_list.append(["EGFR","-","Exon 19 Deletion","Exon 19 Deletion"]+["no"])
                                else:
			    	    print >>snv_out,"\t".join([gene]+snv_sites[gene][site]+["no"]+["-"]*len(need_data))	
				    snv_out_list.append([gene]+snv_sites[gene][site]+["no"])
	snv_out.close()
	return snv_out_list,snv_detacted

def read_fusion_txt(fusion_txt,min_su):
	d=cl.defaultdict(list)
	fusion = open(fusion_txt)
	header = fusion.readline().strip("#\n").split("\t")
	data_d=cl.namedtuple("fusion",header)
	r=lambda x:re.search(",(\S+);",x).group(1)
        for x in fusion:
                x=x.strip().split("\t")
                if len(x) < len(header): continue
                x_d=data_d(*x)
		m=re.findall("(\S+)\s+([^(]+)\([^)]+\)-(\S+)\s+([^(]+)\(",x_d.fusion_format)
		try:
			su,sr,pe=x_d.sample_support.split(",")[0].split(":")
			sr_r=r(x_d.ratio_sr)
			pe_r=r(x_d.ratio_pe)
		except:
			su,sr,pe="000"
			sr_r,pe_r="--"
		if m  and int(su) >min_su and int(sr)*int(pe) >0:
			su_r = sum(map(float,[sr_r,pe_r]))
			trans = lambda x:"{:.3f}".format(float(x)*100)
			d[m[0]] = ["yes",su,sr,pe,trans(su_r),trans(sr_r),trans(pe_r)]

	return d

def make_fusion_results(fusion_txt,fusion_sites,min_su,sample_code,results_dir):
#	fusion_txt = glob.glob(os.path.join(results_dir,"04.fusion_detect","*_lumpy_fusion.desision.vcf"))
	fusion_report = open(os.path.join(results_dir,sample_code + "_fusion_report.xls"),"w")
	if not fusion_txt:
		print >>fusion_report, "can not find fusion_out!"
		fusion_report.close()
		return "can not find fusion_out",0
	d = read_fusion_txt(fusion_txt,min_su)
	fusion_out_list = []
	fusion_detected = cl.defaultdict(int)
	print >>fusion_report,"\t".join(["gene","fusion","pattern","detected","su","sr","pe","ratio_su","ratio_sr","ratio_pe"])
	target_fusion = sorted(fusion_sites.keys())
	detected_fusion = d.keys()  
	for fusion in target_fusion:
		for df in detected_fusion : #df == detected_fusion
			if "-".join([df[0],df[2]]) == fusion:
				if "EML4-ALK"== fusion and df[3][2:] in {"19","20"}:
					pattern = "E{0}A{1}".format(df[1][2:],df[3][2:])
					print >>fusion_report,"{0}\t{1}\t{2}\t{3}".format(fusion_sites[fusion],fusion,pattern,"\t".join(d[df]))
					fusion_out_list.append([fusion_sites[fusion],fusion,pattern,"-"]+["yes"])
					fusion_detected[fusion] +=1
				else:
					pattern = "{0}-{1}-{2}-{3}".format(df[0],df[1][2:],df[2],df[1],df[3][2:])
					print >>fusion_report,"{0}\t{1}\t{2}\t{3}".format(fusion_sites[fusion],fusion,pattern,"\t".join(d[df]))
					fusion_out_list.append([fusion_sites[fusion],fusion,pattern,"-"]+["yes"])
					fusion_detected[fusion] +=1
		if fusion_detected[fusion] < 1:
			print >>fusion_report,"{0}\t{1}\t{2}\t{3}\t{4}".format(fusion_sites[fusion],fusion,"-","no","\t".join(["-"]*6))
			fusion_out_list.append([fusion_sites[fusion],fusion,"-","-"]+["no"])
	all_detected_fusion = sum(fusion_detected.values())
	return fusion_out_list,all_detected_fusion
				
def get_parser():
	parser= argparse.ArgumentParser(description="collect report data.")
	add=parser.add_argument
	add('-d',metavar='results_dir',default = os.getcwd(), help = 'path to results dir')
	add('-f',metavar='fusion_txt',help = 'fusion_txt')
	add('-r',metavar='report',help = 'report.xls')
	add('-p',metavar='project_dir',help = 'project_dir.')
	add('-c',metavar='config_file.', default = "/gpfs/home/shizhenyu/project/20180503_baozheng_data/qc_report/180723/results/temp.yaml",help='config_file.')
	add('-s',metavar='sample_code',type = str,default = "unnamed",help = 'sample_code')
	args=parser.parse_args()
	return args

	
if __name__ == "__main__":
	args = get_parser()
	config_d = yaml.load(open(args.c))
	all_report_d = {"qcFailed":0,"snvDetected":0,"fusionDetected":0,"fusionList":[],"snvList":[]}	
	qc_fail_list = make_qc_report(args.p,config_d["qc_cutoff"],config_d["vip_sites"],args.s,args.d)
	if qc_fail_list :
		all_report_d["qcFailed"] = len(qc_fail_list)
	all_report_d["fusionList"],all_report_d["fusionDetected"] = make_fusion_results(args.f,config_d["fusion_sites"],config_d["min_su"],args.s,args.d)
	all_report_d["snvList"],all_report_d["snvDetected"] = make_snv_report(args.r,config_d["vip_sites"],args.s,args.d)
	out_json = open(os.path.join(args.d, args.s + ".report.json"),"w")	
	json.dump(all_report_d,out_json)
