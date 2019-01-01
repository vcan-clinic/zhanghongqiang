#coding=utf-8
import sys
import re
import collections 
import gzip
from indel_anno_3 import indel_anno
import yaml
import os
class annotation():
	def __init__(self,reformed_vcf,multianno_txt,config_dict):
		self.indel_anno=indel_anno(config_dict["cds_fasta"])
		self.vip_sites=config_dict["vip_sites"]
		self.filter_args=config_dict["filter_args"]
		self.reformed_vcf=self.read_file(reformed_vcf)
		self.multianno_txt=self.read_file(multianno_txt)
		self.name=multianno_txt.split("/")[-1]
		self.position_anno=collections.defaultdict(dict)
		self.vcf_filter()
		self._get_position_anno()

	def tranvar_anno(self,NM,c_anno):
		p_extract = lambda x:x.split("\t")[4].split("/")[2] #
		m = re.match(r'(c\.)?([\d*+-]+)(_([\d*+-]+))?(\.)?(del([atgcnATGCN\d]*))?(ins([atgcnATGCN]*))?(([atgcnATGCN?]*)>([atgcnATGCN?]*))?(dup([atgcnATGCN\d]*))?$', c_anno)
		if  not m:
			m1 = re.match("c\.(\d+)_(\d+)([atgcnATGCN]+)$",c_anno)
			if m1 : 
				c_anno = "c.{}_{}delins{}".format(*m1.groups())		
		if NM.startswith("ENST"):
			out =  os.popen("transvar canno -i '{NM}:{c_anno}' --ensembl --noheader".format(**locals())).read()
		else:
			out =  os.popen("transvar canno -i '{NM}:{c_anno}' --ucsc --noheader".format(**locals())).read()
		try:
			p_anno =  out.split("\t")[4].split("/")[2]
		except:
			p_anno = ""
		if p_anno in {"","."}:
			return ""
		else:
			return p_anno
			
	def read_file(self,f):
		f_open=gzip.open(f) if f.endswith(".gz") else open(f)
		header=f_open.readline().strip().split("\t")
		for line in f_open:
			line=line.strip().split("\t")
			line_header=header+["info_"+str(x) for x in range(len(line)-len(header))]
			yield dict(zip(line_header,line))
		f_open.close()

	def get_cosmic_id(self,cosmic_data):
		try:
			cosmic_id=re.search("ID=(\S+?);",cosmic_data).group(1)
			return cosmic_id
		except :
			return cosmic_data

	def get_tp_anno(self,line):
		tp_anno=line["AAChange.refGene"]
		all_tp=[x.split(":") for x in tp_anno.split(",")] #[[gene,tp,exon,c_anno,(p_anno)],[]]
		if len(all_tp[0])>=4:
			c_anno=[x[3] for x in all_tp]
			tp=[x[1] for x in all_tp]
			exon=[x[2] for x in all_tp]	
			if len(line["Ref"]) >1 or len(line["Alt"]) >1:
				p_anno = []
				for x in all_tp:
					p_transvar = self.tranvar_anno(x[1],x[3])
					if p_transvar:
						p_anno.append(p_transvar)
					else:
						p_anno.append(self.indel_anno.get_posi(x[3],x[1],line["Ref"]))
			else:
				p_anno=[x[4] if len(x)>=5 else "" for x in all_tp]
		else:
			c_anno,p_anno,tp,exon=[],[],[],[]
		return c_anno,p_anno,tp,exon

	def str_list(self,l):
		if type(l)!= str:
			try:
				return ",".join(l)
			except:
				return ""
		else:
			return l 	
	
	def cal_tag(self,line,mut_id):
		if line["Gene.refGene"] in self.vip_sites:
			if mut_id in self.vip_sites[line["Gene.refGene"]]: #是否目标位点
				return "vip"
			elif line["Gene.refGene"] == "EGFR" and  mut_id in self.vip_sites["EGFR"]["19del"]: #是否在cosmic的19del
				return "vip_19del"
			elif mut_id[0] == "chr7" and  (55242415 <= int(line["Start"]) <= 55242513 or 55242415 <= int(line["End"]) <= 55242513) and len(mut_id[2]) > len(mut_id[3]) :#[chr7, 55242465, 55242479, COSM6223, p.E746_A750delELREA]
				return "unknown_19del"
			else:
				return "others"
		else:
			return "others"
			
	def  _get_position_anno(self):
		key=["Chr","info_3","info_5","info_6"]
		for line in self.multianno_txt:
			mut_id=tuple([line[x] for x in key])
			if mut_id in self.position_anno :
				line["cosmic70"]=self.get_cosmic_id(line["cosmic70"])
				line["CDS_Mutation"],line["AA_Mutation"],line["transcripts"],line["exon"]=self.get_tp_anno(line)
				line["tag"]=self.cal_tag(line,mut_id)
				self.position_anno[mut_id]["anno_data"]=line
	
	def filter_rules(self,line,filter_type): #filter_type= "vcf","anno"
		d=self.filter_args
		coverage=lambda x : int(x["alt_depth"]) >d["alt_depth"] and int(x["all_depth"]) >d["all_depth"]
		vf=lambda x : float(line["tumor_ratio"]) >d["tumor_ratio"]# and float(line["tumor_ratio"]) >20*float(line["normal_ratio"])
		db_freq = 0 if line.get("1000g2015aug_all",0) == "." else float(line.get("1000g2015aug_all",0))
		genome_1000=lambda x :db_freq < d["db_freq"]
		def bias(line):
			m=re.search("\d+:(\d+),(\d+)",line["tumor_alt_forward_reverse"])
			if m :
				r1,r2=map(float,m.groups())
			else:
				print line["tumor_alt_forward_reverse"]
				r1,r2=5,5
			return r1 !=0 and r2 !=0 and 0.1 < r1/r2 <10  
		if filter_type == "vcf":
			rules=[coverage,bias,vf]
		elif filter_type == "anno":
			rules=[genome_1000]
		if False in {r(line) for r in rules}:
			return False
		else:
			return True

	def vcf_filter(self):
		key=["chr","pos","ref","alt"]
		#software=["varscan-pair","vardictjava-pair","mutect2-pair"]
		for line in self.reformed_vcf:
			line["normal_ratio"] = "0" if line["normal_ratio"]=="." else line["normal_ratio"]
			mut_id=tuple([line[x] for x in key])
			filter_tag = self.filter_rules(line,"vcf")
			if not filter_tag:
				continue
			if "data" in  self.position_anno[mut_id]:
				self.position_anno[mut_id]["data"].append(line)
			else:
				self.position_anno[mut_id]["data"]=[line]
	
	def print_vcf(self, out_file):
		#out=open(self.name+"_filtered.xls","w")
		out=open(out_file, "w")
	#	temp_d={for x in  }
		software=["varscan-pair","vardictjava-pair","mutect2-pair"]
	#	data_header=['chr', 'pos', 'gene_symbol', 'ref', 'alt', 'homopolymer', 'tumor_ratio', 'normal_ratio', 'genotype', 'alt_depth', 'all_depth', 'control_alt_depth', 'control_all_depth', 'status', 'tumor_alt_forward_reverse', 'normal_alt_forward_reverse', 'decision', '.' ]
	#	data_header_print=['chr', 'pos', 'gene_symbol', 'ref', 'alt', 'homopolymer', 'tumor_ratio', 'normal_ratio', 'genotype', 'alt_depth', 'all_depth', 'control_alt_depth', 'control_all_depth', 'status', 'tumor_alt_forward_reverse', 'normal_alt_forward_reverse', 'decision', 'software' ]
	#	anno_header=["tag",'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'cosmic70', 'avsnp147', 'CLINSIG', 'CLNDBN',"CDS_Mutation","AA_Mutation","exon","transcripts","1000g2015aug_all"]
		data_header=['chr', 'pos', 'ref', 'alt', 'tumor_ratio', 'normal_ratio', 'genotype', 'alt_depth', 'all_depth', 'control_alt_depth', 'control_all_depth', 'status', 'tumor_alt_forward_reverse', 'normal_alt_forward_reverse', 'decision']
		data_header_print=data_header
		anno_header=["Gene.refGene","tag", 'cosmic70',"CDS_Mutation","AA_Mutation","exon","transcripts"]
		print >>out,"\t".join(anno_header+data_header_print)
		for mut_id in self.position_anno:
			if "data" in self.position_anno[mut_id]:
				for line in self.position_anno[mut_id]["data"]:
					data=[line[x] for x in data_header]	
					filter_tag = self.filter_rules(self.position_anno[mut_id].get("anno_data",{}),"anno")
					if not filter_tag :
						continue
			#		if int(line["alt_depth"]) <8:print line 
					anno=[self.str_list(self.position_anno[mut_id].get("anno_data",{}).get(x,"NA")) for x in anno_header]
					if anno[0]=="ERBB2" :anno[0] = "HER2"				
					print >>out, "\t".join(anno+data)
		out.close()
if __name__ == "__main__":
	config_dict=yaml.load(open(sys.argv[3])) #yaml
	a = annotation(sys.argv[1],sys.argv[2],config_dict) #reformed_vcf,multianno_txt
	a.print_vcf(sys.argv[4])#out file name 
