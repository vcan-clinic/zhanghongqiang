import sys
import re
import glob
import os
import gzip
import yaml
import numpy as np
from  collections import defaultdict

class QC_report:

	def __init__(self,path,unit,vipsites):
		self.sample_path=path.rstrip("/")
		self.data=defaultdict(dict)
		self.vipsites=vipsites
		self.site_data=defaultdict(dict)
	def read_Basic_Statistics(self,f):
		d={}
		with open(f) as fo:
			fo.readline()
			for line in fo:
				line=line.strip().split("\t")
				if line:
					d[line[0].strip()]=[re.findall("(\d+)\s*\(\s*([.\d]+)\%\)",x) for x in line[1:]]
		return d
				
	def get_q30(self,sample_type="tumor"):
		trim_Basic_Statistics_txts=glob.glob(os.path.join(self.sample_path,"qc","trim_adapter",sample_type,"*","Basic_Statistics_of_Sequencing_Quality.txt"))
		clean_Basic_Statistics_txts=glob.glob(os.path.join(self.sample_path,"qc","clean",sample_type,"*","Basic_Statistics_of_Sequencing_Quality.txt"))
		raw_q30=0
		clean_q30=0
		total_bases=0
		clean_bases=0 
		p=lambda x:"{:.2f}".format(x*100)
		for txt in  trim_Basic_Statistics_txts:
			txt_data=self.read_Basic_Statistics(txt)
			raw_q30+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][0][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][2][0][0]]))
#			clean_q30+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][1][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][3][0][0]]))
			total_bases+=sum(map(int,[txt_data["Total number of bases"][0][0][0],txt_data["Total number of bases"][2][0][0]]))
#			clean_bases+=sum(map(int,[txt_data["Total number of bases"][1][0][0],txt_data["Total number of bases"][3][0][0]]))
		for txt in  clean_Basic_Statistics_txts:
			txt_data=self.read_Basic_Statistics(txt)
			if not trim_Basic_Statistics_txts:
				raw_q30+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][0][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][2][0][0]]))
				total_bases+=sum(map(int,[txt_data["Total number of bases"][0][0][0],txt_data["Total number of bases"][2][0][0]]))
			clean_q30+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][1][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][3][0][0]]))	
			clean_bases+=sum(map(int,[txt_data["Total number of bases"][1][0][0],txt_data["Total number of bases"][3][0][0]]))
		self.data[sample_type]["Raw_Q30 %"] = p(raw_q30/(total_bases+1))
		self.data[sample_type]["Clean_Q30 %"] = p(clean_q30/(clean_bases+1))
		self.data[sample_type]["Total_bases"] = total_bases
		self.data[sample_type]["clean_bases"] = clean_bases
		if raw_q30 == 0:
			print self.sample_path
	def get_maping_data(self,sample_type="tumor"): 
		dup_dir=glob.glob(os.path.join(self.sample_path,"qc_report","*","dupQC"))
		dedup_dir=glob.glob(os.path.join(self.sample_path,"qc_report","*","dedupQC"))
		coverage_report_dup=self.read_coverage_report(os.path.join(dup_dir[0],"coverage.report")) if dup_dir else defaultdict(str)
		coverage_report_dedup=self.read_coverage_report(os.path.join(dedup_dir[0],"coverage.report")) if dedup_dir else defaultdict(str)
		depth_gz=os.path.join(dedup_dir[0],"depth.tsv.gz") if dedup_dir else ""
		self.data[sample_type]["depth_shift"],self.data[sample_type]["Median_depth"]=self.cal_depth_shift(depth_gz)
		self.site_data[sample_type]=self.depth_stats(depth_gz)
		self.data[sample_type]['Length of region'] = coverage_report_dup["[Target] Len of region"]
		self.data[sample_type]['Fraction of mapped bases(%)'] = coverage_report_dup['[Total] Fraction of Mapped Data(Mb)']
		self.data[sample_type]['Fraction of PCR duplicate(%)'] = coverage_report_dup['[Total] Fraction of PCR duplicate reads']
		self.data[sample_type]['Fraction of target PCR duplicate(%)'] =\
		    round(  (float(coverage_report_dup["[Target] Target Data(Mb)"]) - float(coverage_report_dup['[Target] Target Data Rmdup(Mb)']))/float(coverage_report_dup["[Target] Target Data(Mb)"])*100, 2)
		self.data[sample_type]['Capture efficiency(rmdup,%)'] = round(float(coverage_report_dup['[Target] Target Data Rmdup(Mb)'])/float(coverage_report_dup['[Total] Mapped Data(Mb)']) ,4)*100
		self.data[sample_type]['Capture efficiency(%)'] = coverage_report_dup['[Target] Fraction of Target Data in mapped data']
		self.data[sample_type]['Average depth'] = coverage_report_dup['[Target] Average depth']
		self.data[sample_type]['Average depth(rmdup)'] = coverage_report_dup['[Target] Average depth(rmdup)']
		self.data[sample_type]['Coverage (>0x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>0x)']
		self.data[sample_type]['Coverage (>=30x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>=30x)']
		self.data[sample_type]['Coverage (>=500x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>=500x)']
		self.data[sample_type]['Coverage (>=1000x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>=1000x)']
		self.data[sample_type]['Coverage (>=2000x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>=2000x)']
		self.data[sample_type]['Fraction of properly paired(%)'] = coverage_report_dup['[Total] Fraction of Properly paired']
		self.data[sample_type]['Forward strand reads'] = coverage_report_dup['[Total] forward strand reads']
		self.data[sample_type]['Backward strand reads'] = coverage_report_dup['[Total] backward strand reads']
		read1_trim_rate,read2_trim_rate = self.parseTrimResult(sample_type)
		self.data[sample_type]["Fq1 trim adapter ratio(%)"] = read1_trim_rate
		self.data[sample_type]["Fq2 trim adapter ratio(%)"] = read2_trim_rate
		self.data[sample_type]["Insert Size"] = self.get_insert_size(os.path.join(dedup_dir[0],"insertsize.plot"))
		self.data[sample_type]["hot_spot_min_depth"]=min(self.site_data[sample_type].values()) if self.site_data[sample_type] else "NA"
		PM=re.findall("(RD[^/_]+|PM[^/_]+|WX[^/_]+)",dedup_dir[0])
		cover_rate=self.cover_rate(os.path.join(dedup_dir[0],"depth_distribution.plot"),500)
		header=sorted(self.data[sample_type].keys())
		for i in header:
			j = re.sub("[^\w]+","_",i)
			print "{0}\t{1}\t{2}".format(PM[0],j,self.data[sample_type][i])
		return self.data[sample_type],self.site_data[sample_type] 
	def cal_depth_shift(self,depth_gz):
		if not os.path.isfile(depth_gz): return "NA"
		with gzip.open(depth_gz) as dp:
			dp.readline()
			depth=[x.strip().split("\t") for x in dp]
		depth=map(float,zip(*depth)[2])
		mean=np.mean(depth)
		median=np.median(depth)
		dp_shift=abs(mean-median)/mean*100
		return dp_shift,median

	def get_insert_size(self,insert_plot_file):
		with open(insert_plot_file) as insert_file:
			data = [x.strip().split("\t") for x in insert_file]
		insert_size = max(data,key=lambda x:float(x[2]))[0]
		return insert_size
		 
	
	def parseTrimResult(self,sample_type="tumor"):
	        read1_trim_result = glob.glob(os.path.join(self.sample_path, "qc/trim_adapter",sample_type, "*","*_R1_*_trimming_report.txt"))
	        read2_trim_result = glob.glob(os.path.join(self.sample_path, "qc/trim_adapter",sample_type,"*","*_R2_*_trimming_report.txt"))
	        def _tmp(files):
			total_reads = 0
			reads_with_adapters = 0
			for each in files:
				with open(each) as e_o:
					each_data=e_o.read()
				m=re.search("Total reads processed:\s+([\d,]+)\s+Reads with adapters:\s+([\d,]+)",each_data)
				if m:
					each_total_reads = m.group(1).replace(",","")
					each_reads_with_adapters = m.group(2).replace(",","")
				total_reads +=int(each_total_reads)
				reads_with_adapters += int(each_reads_with_adapters)
			return round(reads_with_adapters/float(total_reads)*100,2)
		read1_trim_rate = _tmp(read1_trim_result)
		read2_trim_rate = _tmp(read2_trim_result)
		return read1_trim_rate,read2_trim_rate
		
	def depth_stats(self,depth_gz):
		if not os.path.isfile(depth_gz): return {}
		sites = defaultdict(list)
		for gene in self.vipsites:
			for site in self.vipsites[gene]:
				if site != "19del":
					sites[site[:2]].append(site)
				else:
					sites[("chr7" ,"55242464")].append("19del")
					
		d={}
		with gzip.open(depth_gz) as dp:
			dp.readline()
			for line in dp:
				line=line.strip().split("\t")
				if tuple(line[:2]) in sites:
					for site in sites[tuple(line[:2])]:
						d[site] = int(line[2])

		return d

	def cover_rate(self,depth_dp,cut_off=300): #depth_distribution.plot
		if not os.path.isfile(depth_dp): return "NA"
		with open(depth_dp) as do:
  			for i in do:
				i=i.strip().split("\t")
				if i[0] == str(cut_off):
					return i[-1]
	def read_coverage_report(self,coverage_report):
		try:
			with open(coverage_report) as cr:
				all_data=re.findall('(\[\S+\]\s+.*?)\t(\S+?)[%]*\n',cr.read())
				return dict(all_data)
		except:
			return defaultdict(str)
	
	def read_file(self,f):
		f_open=gzip.open(f) if f.endswith(".gz") else open(f)
		header=f_open.readline().strip().split("\t")
		for line in f_open:
			line=line.strip().split("\t")
			line_header=header+["info_"+str(x) for x in range(len(line)-len(header))]
			yield dict(zip(line_header,line))
		f_open.close()

	def get_mut_rate(self):
		vcf_gz=glob.glob(os.path.join(self.sample_path,"annotation","reformed.*tool.vcf.gz"))
		key=["chr","pos","ref","alt"]	
		cut_r=[x*0.1 for x in range(2,11)+[20,30,40,50,60,100]]
		d=defaultdict(set)
		if vcf_gz :
			#PM=re.search("(PM\d+)",vcf_gz[0])
			for line in self.read_file(vcf_gz[0]):
				mut_id=tuple([line[x] for x in key])
				ratio=float(line["tumor_ratio"])
				for t in cut_r:
					if ratio < t :
						d[t].add(mut_id)	
			rate="\t".join([str(len(d[x])) for x in cut_r])
			for line in self.read_file(last[0]):
				mut_id=tuple([line[x] for x in key])
				ratio=float(line["tumor_ratio"])
				for t in cut_r:
					if ratio < t :
						d[t].add(mut_id)	
			print "{0}\t{1}\t{2}\t{3}\t{4}".format(self.data["Raw_Q30"],self.data["Clean_Q30"],self.data["Total_bases"],self.data["clean_bases"],rate)
		else:
			return 0
			
		
		
if __name__ =="__main__":
	config_d = yaml.load(open("/gpfs/home/shizhenyu/scripts/mk_report/report.yaml"))
	for i in open(sys.argv[1]):
		i = i.split("\t")
		project_dir = i[0]
		try:
			q = QC_report(project_dir,"",config_d['vip_sites'])
			q.get_q30("tumor")
			qc_d,site_depth_d = q.get_maping_data()
		except:
			print "ERROR!!! "+project_dir
