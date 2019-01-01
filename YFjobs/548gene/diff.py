7d6
< import json
13c12
< 	def __init__(self,path,unit_txt,vipsites):
---
> 	def __init__(self,path):
15,16d13
< 		self.info={"tumor":"","tumor_fastq":[],"normal":"","normal_fastq":[]}
< 		self._sample_unit(unit_txt)
18,19c15
< 	#	self.vipsites=yaml.load(open("/gpfs/home/shizhenyu/project/20180222_baozhengruanjian/anno_test/vip_sites.yaml"))["vip_sites"]
< 		self.vipsites=vipsites
---
> 		self.vipsites=yaml.load(open("/gpfs/home/shizhenyu/project/20180222_baozhengruanjian/anno_test/vip_sites.yaml"))["vip_sites"]
21,33d16
< 		self.get_q30()
< 	def _sample_unit(self,unit_txt):
< 		with open(unit_txt) as ut:
< 			data = [i.strip().split() for i in ut]
< 			samples=zip(*data)[1]
< 			self.info["tumor"] = samples[0]
< 			self.info["normal"] = samples[-1] 
< 			for fastq_prefix,sample in data:
< 				if sample == self.info["tumor"] :
< 					self.info["tumor_fastq"].append(fastq_prefix)
< 				if sample == self.info["normal"]:
< 					self.info["normal_fastq"].append(fastq_prefix)
< 
44,51c27,31
< 	def get_q30(self,sample_type="tumor"):	
< 		fastp_jsons = []
< 		clean_Basic_Statistics_txts= []
< 		for fastq_prefix in self.info[sample_type+"_fastq"]:
< 			clean_Basic_Statistics_txts += glob.glob(os.path.join(self.sample_path,"01.qc",fastq_prefix,"Basic_Statistics_of_Sequencing_Quality.txt"))
< 			fastp_jsons += glob.glob(os.path.join(self.sample_path,"01.qc","fastp",fastq_prefix+"*.json"))
< 		raw_q30_bases=0
< 		clean_q30_bases=0
---
> 	def get_q30(self,sample_type="tumor"):
> 		trim_Basic_Statistics_txts=glob.glob(os.path.join(self.sample_path,"qc","trim_adapter",sample_type,"*","Basic_Statistics_of_Sequencing_Quality.txt"))
> 		clean_Basic_Statistics_txts=glob.glob(os.path.join(self.sample_path,"qc","clean",sample_type,"*","Basic_Statistics_of_Sequencing_Quality.txt"))
> 		raw_q30=0
> 		clean_q30=0
55,58c35,40
< 		for json_f in fastp_jsons:
< 			temp_d=json.load(open(json_f))["summary"]["before_filtering"]
< 			total_bases+=temp_d["total_bases"]
< 			raw_q30_bases+=temp_d["q30_bases"]
---
> 		for txt in  trim_Basic_Statistics_txts:
> 			txt_data=self.read_Basic_Statistics(txt)
> 			raw_q30+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][0][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][2][0][0]]))
> #			clean_q30+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][1][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][3][0][0]]))
> 			total_bases+=sum(map(int,[txt_data["Total number of bases"][0][0][0],txt_data["Total number of bases"][2][0][0]]))
> #			clean_bases+=sum(map(int,[txt_data["Total number of bases"][1][0][0],txt_data["Total number of bases"][3][0][0]]))
61c43,46
< 			clean_q30_bases+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][1][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][3][0][0]]))	
---
> 			if not trim_Basic_Statistics_txts:
> 				raw_q30+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][0][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][2][0][0]]))
> 				total_bases+=sum(map(int,[txt_data["Total number of bases"][0][0][0],txt_data["Total number of bases"][2][0][0]]))
> 			clean_q30+=sum(map(float,[txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][1][0][0],txt_data["Number of base calls with quality value of 30 or higher (Q30+) (%)"][3][0][0]]))	
63,65c48,50
< 		self.data[sample_type]["Raw_Q30 %"] = p(float(raw_q30_bases)/(total_bases+1))
< 		self.data[sample_type]["Clean_Q30 %"] = p(clean_q30_bases/(clean_bases+1))
< 		self.data[sample_type]["Total_bases"] = total_bases 
---
> 		self.data[sample_type]["Raw_Q30"] = p(raw_q30/(total_bases+1))
> 		self.data[sample_type]["Clean_Q30"] = p(clean_q30/(clean_bases+1))
> 		self.data[sample_type]["Total_bases"] = total_bases
67,70c52,56
< 
< 	def get_maping_data(self,sample_type="tumor"): #self.tumor,self.normal
< 		dup_dir=glob.glob(os.path.join(self.sample_path,"02.mapping","base_recal.qc_report_dup",self.info["tumor"]))
< 		dedup_dir=glob.glob(os.path.join(self.sample_path,"02.mapping","base_recal.qc_report_dedup",self.info["tumor"]))
---
> 		if raw_q30 == 0:
> 			print self.sample_path
> 	def get_maping_data(self,sample_type="tumor"): 
> 		dup_dir=glob.glob(os.path.join(self.sample_path,"qc_report","*","dupQC"))
> 		dedup_dir=glob.glob(os.path.join(self.sample_path,"qc_report","*","dedupQC"))
74,76c60
< 		self.data[sample_type]["depth_shift"]=self.cal_depth_shift(depth_gz)
< #		print self.data
< 		depth_gz = os.path.join(dedup_dir[0],"depth.tsv.gz")
---
> 		self.data[sample_type]["depth_shift"],self.data[sample_type]["Median_depth"]=self.cal_depth_shift(depth_gz)
78,79c62
< 		self.site_data[sample_type]
< 		self.data[sample_type]['Length of region'] = coverage_report_dedup["[Target] Len of region"]
---
> 		self.data[sample_type]['Length of region'] = coverage_report_dup["[Target] Len of region"]
90,91c73,74
< 		self.data[sample_type]['Coverage (>=100x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>=100x)']
< 	#	self.data[sample_type]['Coverage (>=1000x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>=1000x)']
---
> 		self.data[sample_type]['Coverage (>=500x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>=500x)']
> 		self.data[sample_type]['Coverage (>=1000x,rmdup,%)'] = coverage_report_dedup['[Target] Coverage (>=1000x)']
100c83
< 		self.data[sample_type]["hot_spot_min_depth"]=min(self.site_data[sample_type].values())
---
> 		self.data[sample_type]["vip_min_depth"]=min(self.site_data[sample_type].values()) if self.site_data[sample_type] else "NA"
103,107c86,105
< #		header=sorted(self.data[sample_type].keys())
< #		for i in header:
< #			print "{0}\t{1}\t{2}".format(PM[0],i,self.data[sample_type][i])
< 		return self.data[sample_type],self.site_data[sample_type]
< 
---
> 		header=sorted(self.data[sample_type].keys())
> 		for i in header:
> 			j = re.sub("[^\w]+","_",i)
> 			print "{0}\t{1}\t{2}".format(PM[0],j,self.data[sample_type][i])
> 		'''
> 		for site in self.data[sample_type]["site_depth"]:
> 			for site_detail in self.vipsites["all_sites"][site[:2]]:
> 				print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(PM[0]+"\t"+sample_type,
> 						self.data[sample_type]["Raw_Q30"],
> 						self.data[sample_type]["Average depth(rmdup)"],
> 						self.data[sample_type]["Coverage (>=500x,rmdup,%)"],
> 						cover_rate,
> 						self.data[sample_type]["Total_bases"]/1000000.0,
> 						self.data[sample_type]["Fraction of mapped bases(%)"],
> 						self.data[sample_type]["depth_shift"],
> 						vip_min_depth,
> 						site_detail[7],
> 						site_detail[6],
> 						self.data[sample_type]["site_depth"][site])
> 		'''
117c115
< 		return dp_shift
---
> 		return dp_shift,median
127,128c125,126
< 	        read1_trim_result = glob.glob(os.path.join(self.sample_path, "01.qc",self.info[sample_type]+"*_R1_*_trimming_report.txt"))
< 	        read2_trim_result = glob.glob(os.path.join(self.sample_path, "01.qc",self.info[sample_type]+"*_R2_*_trimming_report.txt"))
---
> 	        read1_trim_result = glob.glob(os.path.join(self.sample_path, "qc/trim_adapter",sample_type, "*","*_R1_*_trimming_report.txt"))
> 	        read2_trim_result = glob.glob(os.path.join(self.sample_path, "qc/trim_adapter",sample_type,"*","*_R2_*_trimming_report.txt"))
148,155d145
< 		sites = defaultdict(list)
< 		for gene in self.vipsites:
< 			for site in self.vipsites[gene]:
< 				if site != "19del":
< 					sites[site[:2]].append(site)
< 				else:
< 					sites[("chr7" ,"55242464")].append("19del")
< 					
161,163c151,163
< 				if tuple(line[:2]) in sites:
< 					for site in sites[tuple(line[:2])]:
< 						d[site] = int(line[2])
---
> 				if line[0] == "chr7" and  55242465< int(line[1]) < 55242479:
> 					if ("chr7","55242465","\\","\\") in  d:
> 						d[("chr7","55242465","\\","\\")].append(int(line[2]))
> 					else:
> 						d[("chr7","55242465","\\","\\")]=[int(line[2])]	
> 				elif tuple(line[:2]) in self.vipsites["all_sites"] and  tuple(line[:2]) != ('chr7', '55242465'):
> 					for site in  self.vipsites["all_sites"][tuple(line[:2])]:
> 						d[tuple(site[:4])] = int(line[2])
> 
> 			d[("chr7","55242465","\\","\\")]=np.mean(d[("chr7","55242465","\\","\\")])
> 	#	if ("chr12", '25398284', "C", "T") in d : 
> 	#		d[("chr12", '25398284', "C", "A")]=d[("chr12", '25398284', "C", "T")]
> 	#	print d.keys()
217,231c217,227
< 	sample_dir,sample_unit=sys.argv[1:3]
< 	q=QC_report(sample_dir,sample_unit)
< 	q.get_q30("tumor")
< 	q.get_maping_data()
< 
< #	d={"wc":"normal","ffpe":"tumor"}
< #	with open(sys.argv[1]) as l:
< #		for i in l:
< #			sample_dir,panel=i.strip().split("\t")
< 	#		sample_dir=os.path.abspath(sample_dir)
< 	#		s_type=re.findall("_(ffpe|wc)",sample_dir)
< 	#		if s_type:
< 	#			sample_type=d[s_type[0]]
< 	#		else:
< 	#			sample_type="tumor"
---
> #	print "\t".join(["rawq30","cleanq30","total_base","clean_base"]+["cutoff_"+str(x*0.1) for x in range(2,11)+[20,30,40,50,60,1000]])
> 	d={"wc":"normal","ffpe":"tumor"}
> 	with open(sys.argv[1]) as l:
> 		for i in l:
> 			sample_dir,panel=i.strip().split("\t")
> 			sample_dir=os.path.abspath(sample_dir)
> 			s_type=re.findall("_(ffpe|wc)",sample_dir)
> 			if s_type:
> 				sample_type=d[s_type[0]]
> 			else:
> 				sample_type="tumor"
233,234c229,230
< 		#	try:
< 	#		if 1:
---
> 			try:
> 		#	if 1:
236,238c232,234
< 	#			q=QC_report(sample_dir)
< 	#			q.get_q30(sample_type)
< 	#			q.get_maping_data(sample_type)
---
> 				q=QC_report(sample_dir)
> 				q.get_q30(sample_type)
> 				q.get_maping_data(sample_type)
241,242c237,238
< 		#	except:
< 	#			pass
---
> 			except:
> 				pass
