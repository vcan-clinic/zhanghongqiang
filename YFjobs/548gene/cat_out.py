import sys,re,glob,os
def read_sample(f):
	try:
		fo = open(f)
	except:
		print f+" can not read!"
		return []
	l = [i.strip().split("\t") for i in fo]
	fo.close()
	return l
	

def cat_out(data_dir,prefix,header_yn = 1):
	f = open(data_dir.split("/")[-1]+"_"+prefix,"w")
	outbox = []
	header = [] 
	for i in glob.glob(os.path.join(data_dir,"*"+prefix)):
		i=i.strip()
		name = re.sub(prefix+"$","",os.path.split(i)[1])
		l = read_sample(i)
		if header_yn:
			header = ["sample_name"]+l[0]
			body = [[name]+x for x in l[1:]]
		else:
			header = ["sample_name"]
			body = [[name]+x for x in l]
		out_print = "\n".join(["\t".join(x) for x in body ])
		outbox.append(out_print)
	print >>f,"\t".join(header)
	print >>f,"\n".join(outbox)
	f.close()

datapath = sys.argv[1]
cat_out(datapath,".report.xls")	
cat_out(datapath,"_fusion_report.xls")	
cat_out(datapath,"_qc_report.xls")	
cat_out(datapath,"_snv_report.xls")	




		

	
