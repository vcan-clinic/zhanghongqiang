#!/usr/bin/python
#encoding: utf-8
import sys,re,os
from Bio import SeqIO
from Bio.Seq import Seq
	
class indel_anno():
	n=0
	def __init__(self,cds_fa=""):
	#	self.data_file=data_file
		if not cds_fa :
			cds_fa="/gpfs/software/pipeline/pipeline/clinical_examination_script/NM_cds.fa"
		self.seq_d=self.get_seq_d(cds_fa)		
	#self.get_posi(self.data_file)

	def read_table(self,f):
	        with open(f) as f_o:
	                l=[i.rstrip("\n").split("\t") for i in f_o]
	        return l
	
	def get_seq_d(self,fasta):
		return {i.id:[i.seq,i.description.split()[1]] for i in SeqIO.parse(fasta,"fasta")}
	
	def get_posi(self,c_anno,NM,ref): #所需数据，CDS_Mutation,NM,ref
		m0=re.match("c.(\d+_\d+)([ATCG]*)(>|del|dup|ins)([ATCG]*)",c_anno)
		m1=re.match("c.(\d+_\d+)([ATCG]+)",c_anno)
		if m0 or m1:
			self.n+=1
			if m0 :
				posi,s1,m_tp,s2=m0.groups()
			else:
				posi,s2=m1.groups()
				m_tp=""
			out=self.trans_p(NM,posi,"multi_m",s2,m_tp,ref)
			if len(out)<9:
				print "indel{0}\t{1},{2}\tannotation failed\t{3}".format(self.n,NM,c_anno,out[0])
				return "annotate_AA_mutation_failed"
			else:
				aa,aa_m,aa_pos_start,aa_pos_end,aa_aj,aa_m_aj,aa_pos_start_aj,aa_pos_end_aj,tag,aa_5,aa_3=out
				try:
					p=self.p_code(*out)
				except:
					print "\t".join(map(str,out))
					p="NA"
				mut_data='{0}(ref_AA:{1},start:{3},end:{4},mut_AA:{2})'.format(p,aa,aa_m,aa_pos_start,aa_pos_end)
				print "indel{0}\t{1},{2}\tannotation succeed\t{3}".format(self.n,NM,c_anno,mut_data)
				return mut_data
		else:
			pass
	 		
	
	def p_code(self,aa,aa_m,aa_pos_start,aa_pos_end,aa_aj,aa_m_aj,aa_start_aj,aa_end_aj,tag,aa_5,aa_3):
		if tag == "fs":
			p="p.{0}{1}fs".format(aa_aj[0],aa_start_aj) if aa_aj else "p.{0}{1}fs".format(aa_3,aa_start_aj)
		elif tag == "=":
			p="p.(=)"
		elif tag == ">": #单氨基酸替换1v1,不加">",如p.A10B;多v1,1v多，或者多V多，都算是indel，写>或者delins 如p.A10>BB或者p.A10_A11>BB
			s_pos=aa_m_aj.find("*")
			aa_m_stop_adj=aa_m_aj[:s_pos+1] if s_pos>=0 else aa_m_aj
		#	if aa_m_aj.find("*") == -1:
		#		left_aa=aa_aj[0]
		#		right_aa=aa_aj[-1]
		#		p="p.{2}{0}_{3}{1}>{4}".format(aa_start_aj,aa_end_aj,left_aa,right_aa,aa_m_aj) if aa_start_aj != aa_end_aj else "p.{2}{0}{3}".format(aa_start_aj,aa_end_aj,aa_aj,aa_m_aj)  
		#	elif len(aa_aj)==1 and len(aa_m_aj)>1:
			if len(aa_aj)==1 and len(aa_m_stop_adj)>1:
				p="p.{1}{0}>{2}".format(aa_start_aj,aa_aj,aa_m_stop_adj)
			else:
				left_aa=aa_aj[0]
				right_aa=aa_aj[-1]
				p="p.{2}{0}_{3}{1}>{4}".format(aa_start_aj,aa_end_aj,left_aa,right_aa,aa_m_stop_adj) if aa_start_aj != aa_end_aj else "p.{1}{0}{2}".format(aa_start_aj,aa_aj,aa_m_stop_adj)
		elif tag == "ins":
			insert_posi=aa_m.rfind(aa_m_aj)-1
			if aa_aj=="":
	                        left_aa=aa_5
	                        right_aa=aa_3
			elif 0 <= insert_posi <=len(aa_aj)-2 : #insert_posi 的位置是ins的前一位,+1是后一位
				left_aa=aa_aj[insert_posi] 
				right_aa=aa_aj[insert_posi+1]   #有些情况是插入到后面位置，比如L 变成了 PAL
			elif insert_posi == -1 :
				left_aa=aa_5
				right_aa=aa_aj[0]
			elif aa_aj=="":
				left_aa=aa_5
				right_aa=aa_3
			else:
				left_aa=aa_aj[insert_posi]
				right_aa=aa_3
			p="p.{2}{0}_{3}{1}ins{4}".format(aa_start_aj,aa_end_aj,left_aa,right_aa,aa_m_aj)
		elif tag == "dup":
			p="p.{2}{0}_{3}{1}dup".format(aa_start_aj,aa_end_aj,aa_m_aj[0],aa_m_aj[-1]) if aa_start_aj != aa_end_aj else "p.{1}{0}dup".format(aa_start_aj,aa_m_aj)  #aa_start_aj != aa_end_aj 意味着dup片段大于1
		elif tag == "del":
			left_aa=aa_m_aj[0]
			right_aa=aa_m_aj[-1]
			p="p.{2}{0}_{3}{1}del{4}".format(aa_start_aj,aa_end_aj,left_aa,right_aa,aa_m_aj) if aa_start_aj != aa_end_aj else "p.{1}{0}del{1}".format(aa_start_aj,aa_m_aj)
		else:
			p="tag not right"
		return p
	
	def clean_flank(self,short_aa,long_aa):
		short_aa_fake_l=short_aa+"#"*(len(long_aa)-len(short_aa))
		# short_aa_fake_r="#"*(len(long_aa)-len(short_aa))+short_aa
		for n,x in enumerate(short_aa_fake_l):
			if long_aa[n] != x :
				break
		for m,x in enumerate(short_aa_fake_l[::-1]):
			if long_aa[-(m+1)]!= x :
				m=-m
				break
	 	return short_aa[n:len(short_aa)+m],long_aa[n:len(long_aa)+m],n,m
	
	def cmp_AA(self,aa,aa_m,tag):
		if tag == "fs":
			if aa == aa_m:
				aa_aj,aa_m_aj,n,m="","",1,0
			elif len(aa_m) > len(aa):
				aa_aj,aa_m_aj,n,m = self.clean_flank(aa,aa_m)
			elif len(aa_m) <= len(aa):
				aa_m_aj,aa_aj,n,m = self.clean_flank(aa_m,aa)
			return "fs",aa_aj,aa_m_aj,n,m
		elif len(aa_m) > len(aa):
			ins=re.findall('([A-Z*]*?)'+'([A-Z*]*)'.join(list(aa))+'([A-Z*]*)',aa_m)
			if ins:
				inserts=[x for x in ins[0] if x]
				if len(inserts) == 1:
					posi=aa_m.rfind(inserts[0])-1  #想要的位置是插入序列的前一位，所以要-1
					if posi+1-len(inserts[0]) >=0 and aa_m[posi+1-len(inserts[0]):posi+1]==inserts[0]:
						return "dup",aa,inserts[0],posi+1-len(inserts[0]),posi
					else:
						return "ins",aa,inserts[0],posi,posi+1
				else:
					aa_aj,aa_m_aj,n,m = self.clean_flank(aa,aa_m)
					return ">",aa_aj,aa_m_aj,n,m
			else:
				aa_aj,aa_m_aj,n,m = self.clean_flank(aa,aa_m)
				return ">",aa_aj,aa_m_aj,n,m
		elif len(aa_m) < len(aa):
			dels=re.findall('([A-Z*]*?)'+'([A-Z*]*)'.join(list(aa_m))+'([A-Z*]*)',aa)	
			if dels:
				deletions=[x for x in dels[0] if x]
				if len(deletions) == 1:	
					posi=aa.rfind(deletions[0]) 
					return "del",aa,deletions[0],posi,posi-1+len(deletions[0])
				else:
					aa_m_aj,aa_aj,n,m = self.clean_flank(aa_m,aa)
					return ">",aa_aj,aa_m_aj,n,m
			else:
				aa_aj,aa_m_aj,n,m = self.clean_flank(aa,aa_m)
				return ">",aa_aj,aa_m_aj,n,m
		elif aa_m == aa:
			return "=",aa,aa_m,0,0
		else:
			aa_aj,aa_m_aj,n,m = self.clean_flank(aa,aa_m)
			return ">",aa_aj,aa_m_aj,n,m
				 
	def trans_p(self,t_id,posi,tp,s2,m_tp,ref):
		if tp == "ins":
			start,end=map(int,posi.split("_"))
			codon_1=(start-1)//3*3+1
			seq=self.seq_d[t_id][start-1:end]
			aa=self.seq_d[t_id][codon_1-1:codon_1+2].translate()
			return seq,aa
		elif tp == "multi_m":
			start,end=map(int,posi.split("_"))
			if m_tp=="ins":
				end-=1
			codon_1=(start-1)//3*3+1 
			codon_2=(end-1)//3*3+1
			
			if t_id not in self.seq_d or len(self.seq_d[t_id][0])< codon_2 :	
				return ("missing transcript seq.",)
			if self.seq_d[t_id][1]=="-":
				c_ref=self.seq_d[t_id][0][start-1:end].reverse_complement()
			else:
				c_ref=self.seq_d[t_id][0][start-1:end]
			if (m_tp not in ["del","ins"]) and (c_ref not in [ref[:len(c_ref)],ref[-len(c_ref):] ]):
				return ("transcript seq not right.",)
			seq=self.seq_d[t_id][0][start-1:end].upper()
			seq_code=self.seq_d[t_id][0][codon_1-1:codon_2+2]
			seq_m=self.seq_d[t_id][0][codon_1-1:start-1]+s2 if m_tp !="ins" else self.seq_d[t_id][0][codon_1-1:start]+s2 #左补全mutation的密码子
			seq_m_code=seq_m+self.seq_d[t_id][0][end:end+{1:2,2:1,0:0}[len(seq_m)%3]] if  len(self.seq_d[t_id][0]) >= end+{1:2,2:1,0:0}[len(seq_m)%3] else seq_m[:3]  #右补全mutation的密码子 
			if len(seq_code)%3 !=0 or len(seq_m_code)%3 !=0:
				return ("mutation on the transcript end.",)
			aa=seq_code.translate()   #原AA序列
			aa_m=seq_m_code.translate() #突变AA序列
			try:
				aa_5=self.seq_d[t_id][0][codon_1-4:codon_1-1].translate() #原AA序列上一位AA
				aa_3=self.seq_d[t_id][0][codon_2+2:codon_2+5].translate() #原AA序列下一位AA
			except:
				aa_5="?"
				aa_3="?"
			aa_pos_start=(start-1)//3+1 
			aa_pos_end=(end-1)//3+1
			if m_tp !="ins" and (end+1-start-len(s2))%3 !=0:
				tag="fs"
			elif m_tp=="ins" and len(s2)%3 !=0:
				tag="fs"
			else:
				tag=""
			tag,aa_aj,aa_m_aj,n,m=self.cmp_AA(str(aa),str(aa_m),tag) #根据突变情况确定是del 还是ins 、普通突变、移码突变
			if tag in {"fs",">"} :
				aa_pos_start_aj=aa_pos_start+ n
				aa_pos_end_aj=aa_pos_end+m
			else:
				aa_pos_start_aj=aa_pos_start+n
				aa_pos_end_aj=aa_pos_start+m
			return aa,aa_m,aa_pos_start,aa_pos_end,aa_aj,aa_m_aj,aa_pos_start_aj,aa_pos_end_aj,tag,aa_5,aa_3
		else:
			return ()
		
if  __name__=="__main__":
	a=indel_anno()	
	header=["Samplecode","tag",'Class','Reported','Score', 'GENE532', 'Software', 'GENE', 'CDS_Mutation', 'AA_Mutation', 'Exon', 'tumor_ratio', 'normal_ratio', 'Leve', 'QC', 'Likely_germline', 'NM', 'status', 'target_chemo_mark', 'COSMIC', 'COSMIC_somatic_status', 'ExonicFunc', 'alt_depth', 'all_depth', 'control_alt_depth', 'control_all_depth', 'tumor_alt_forward_reverse', 'normal_alt_forward_reverse', 'homopolymer', 'dbsnp', 'ALL_1000G', 'EAS_1000G', 'ExAC_Freq', 'ExAC_EAS', 'ESP6500si_ALL','ClinVar_SIG', 'ClinVar_DIS','Pathogenicity','soft_count', 'Freebayes_num', 'Vardict_num','Mutation_db', 'DISEASE', 'FDACFDA1', 'FDACFDA2', 'DRUGS', 'INFORMATION','COSMIC_ID', 'chr', 'pos', 'gene_symbol', 'ref', 'alt', 'SIFT_pred', 'MutationAssessor_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'FATHMM_pred', 'PROVEAN_pred', 'fathmm_MKL_coding_pred', 'MetaSVM_pred', 'MetaLR_pred','info1','info2','info3','info4','info5','info6','info7','info8','info9','info10','Chr','Position','Ref','A','T','C','G',
'N']
	raw_data=open(sys.argv[1])
 #      new_file_name=raw_data_file.replace('.xls','.xls.temp')
 #      out_file = open(new_file_name,"w")
#raw_data_header=raw_data.readline()
#print >>out_file,raw_data_header.strip("\n")
	n=0
	for j in raw_data:
		j=j.strip("\n").split("\t")
		i=j+[""]*(len(header)-len(j))
		data=dict(zip(header,i))
		a.get_posi(data["CDS_Mutation"],data["NM"],data["ref"])
			
