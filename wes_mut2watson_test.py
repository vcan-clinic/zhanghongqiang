import argparse
import os
import gzip
import re
import glob
import time
def header():
    vcf_head = '''##fileformat=VCFv4.1
##INFO=<ID=WfG_variant_origin,Number=1,Type=String,Description="variant origin,value is somatic or germline">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the event in the tumor">
##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">'''
    vcf_head_contig='''##contig=<ID=chr1,assembly=hg19,length=249250621>
##contig=<ID=chr10,assembly=hg19,length=135534747>
##contig=<ID=chr11,assembly=hg19,length=135006516>
##contig=<ID=chr12,assembly=hg19,length=133851895>
##contig=<ID=chr13,assembly=hg19,length=115169878>
##contig=<ID=chr14,assembly=hg19,length=107349540>
##contig=<ID=chr15,assembly=hg19,length=102531392>
##contig=<ID=chr16,assembly=hg19,length=90354753>
##contig=<ID=chr17,assembly=hg19,length=81195210>
##contig=<ID=chr18,assembly=hg19,length=78077248>
##contig=<ID=chr19,assembly=hg19,length=59128983>
##contig=<ID=chr2,assembly=hg19,length=243199373>
##contig=<ID=chr20,assembly=hg19,length=63025520>
##contig=<ID=chr21,assembly=hg19,length=48129895>
##contig=<ID=chr22,assembly=hg19,length=51304566>
##contig=<ID=chr3,assembly=hg19,length=198022430>
##contig=<ID=chr4,assembly=hg19,length=191154276>
##contig=<ID=chr5,assembly=hg19,length=180915260>
##contig=<ID=chr6,assembly=hg19,length=171115067>
##contig=<ID=chr7,assembly=hg19,length=159138663>
##contig=<ID=chr8,assembly=hg19,length=146364022>
##contig=<ID=chr9,assembly=hg19,length=141213431>
##contig=<ID=chrX,assembly=hg19,length=155270560>
##contig=<ID=chrY,assembly=hg19,length=59373566>
##contig=<ID=chrM,assembly=hg19,length=16569>'''
    return (vcf_head,vcf_head_contig)

def poly_data(name,fre):
    vv='fail'
    if fre=='.' or fre =='0' or float(fre) <=0.01:
        vv='pass'
    return vv

def in_cosmic_times(infor):
    '''ID=COSM1140131,COSM553;OCCURENCE=15(lung),2(thyroid),2(haematopoietic_and_lymphoid_tissue),1(stomach),1(small_intestine),1(breast),4(skin),3(biliary_tract),52(large_intestine),1(kidney),2(pancreas),2(urinary_tract),1(soft_tissue),3(endometrium),1(salivary_gland)'''
    if infor==".":
        tag='fail'
    else:
        infor=re.split(';',infor)[1]
        times= re.findall('(\d+)'+'\(',infor)
        total=sum(([int(x) for x in times]))
        tag='fail'
        if total >=4:
            tag='pass'
    return tag

def QC_tissue(pair,infos=dict,onlymutect2=True):
    if pair==True:
        ##add germline
        dp_cut=20
        rt_somatic=1
        dic={}
        qc_rec={}
        info=infos
        # print info
        #GT:DP:VD:ALD
        for i in info:
            # print i
            if onlymutect2: ##only mutect2 considered
                if i=='Software':
                    qc = 'pass' if info['Software'].find('mutect2') >=0 else 'fail'
                    qc_rec[qc]=1
            else:
                if i=='soft_count':###still use mutect2 as qc's software
                    if int(info[i]) >= 2:
                        qc = 'pass'
                        #     qc='pass' if info['Software'].find('mutect2') >=0 else 'fail'
                    else:
                        qc = 'fail'
                    qc_rec[qc] = 1
            if i=='all_depth' or i== 'control_all_depth':
                # if info[i] == "." : info[i] = 0
                # print int(float(info[i]))
                qc='pass' if int(float(info[i]))>=dp_cut else 'fail'
                dp_germ='pass' if int(float(info[i]))>=dp_cut else 'fail'
                qc_rec[qc] = 1
                #print 2,qc
            if i=='tumor_ratio%':
                print 1111,info[i]
                # if info[i] == ".": info[i] = 0
                qc='pass' if (float(info[i]) >=rt_somatic) or\
                             (float(info[i]) >= 10 and
                              (int(re.split("[:,]",info["tumor_alt_forward_reverse"])[0])>=2 and
                               int(re.split("[:,]",info["tumor_alt_forward_reverse"])[1])>=1 and
                               int(re.split("[:,]",info["tumor_alt_forward_reverse"])[2])>=1)
                              ) else 'fail'
                # qc="pass" if float(info[i]) >= 10 and  (int(re.split("[:,]",info["tumor_alt_forward_reverse"])[0])>=2 and int(re.split("[:,]",info["tumor_alt_forward_reverse"])[1])>=1 and int(re.split("[:,]",info["tumor_alt_forward_reverse"])[2])>=1) else "fail"
                qc_rec[qc] = 1
                #print 3,qc
            if i=='normal_ratio%':
                print i
                print 2222,info["tumor_ratio%"]
                print re.split("[:,]", info["tumor_alt_forward_reverse"])
                # if info[i] == ".": info[i] = 0
                # qc = 'pass' if float(info[i]) <= 1 else 'fail'
                if (
                        float(info[i])<=1 and
                        float(info[i])==0 and
                        float(info["tumor_ratio%"]) >=1 and
                        (int(re.split("[:,]",info["tumor_alt_forward_reverse"])[0])>=4and
                         int(re.split("[:,]",info["tumor_alt_forward_reverse"])[1])>=2 and
                         int(re.split("[:,]",info["tumor_alt_forward_reverse"])[2])>=2))or(
                        (float(info[i])<=1 and float(info[i])!=0 and float(info["tumor_ratio%"]) >=1 and
                         (float(info["tumor_ratio%"])/float(info[i])>=10) and
                         (int(re.split("[:,]",info["tumor_alt_forward_reverse"])[0])>=6and
                          int(re.split("[:,]",info["tumor_alt_forward_reverse"])[1])>=3 and
                          int(re.split("[:,]",info["tumor_alt_forward_reverse"])[2])>=3))):
                    qc="pass"
                else:
                    qc="fail"

                # qc_rec[qc] = 1

                # if float(info[i])!=0 and
                #  float(info["tumor_ratio%"]) >=1 and
                #  float(float(info[i])/float(info["tumor_ratio%"])>=10)  and
                #  (int(re.split("[:,]",info["tumor_alt_forward_reverse"])[0])>=6and
                #  int(re.split("[:,]",info["tumor_alt_forward_reverse"])[1])>=3 and
                #  int(re.split("[:,]",info["tumor_alt_forward_reverse"])[2])>=3):
                #     qc="pass"
                # else:
                #     qc="fail"
                qc_rec[qc] = 1
                # print qc_rec[qc]
                #print 4,qc
            # if i=='tumor_alt_forward_reverse':
            #     rd=re.split('[:,]',info[i])
            #     # print rd
            #     qc='pass' if (int(rd[1])>0 and int(rd[2])>0) else 'fail'
            #     qc_rec[qc] = 1
                #print 5,qc
            if i in ['ExAC_Freq','1000G_ALL','ESP6500si_ALL','1000G_EAS','ExAC_EAS']:
                tag=poly_data(i,info[i])
                qc_rec[tag]=1

        qc_germ='fail'
        comb='het'
        if info['Software'].find('vardictjava') >= 0:
            if int(float(info['all_depth']) )>=20 and int(float(info['control_all_depth'])) >=20 and float(info['normal_ratio%'])>=1:
                rt_t=float(info['tumor_ratio%'])
                rt_s=float(info['normal_ratio%'])
                poly='pass'
                for i in ['ExAC_Freq','1000G_ALL','ESP6500si_ALL','1000G_EAS','ExAC_EAS']:
                    tag=poly_data(i,info[i])
                    # print i,info[i],tag
                    if tag=='fail':
                        poly='fail'
                        break
                cos=in_cosmic_times(info['COSMIC'])
                #print poly
                #print cos,rt_t,rt_s
                #################################
                if rt_t >= 20 and rt_s >= 20  :#and rt_t/rt_s>=0.25 and rt_t/rt_s<=4:
                    if poly=='pass' or (poly=='fail' and cos=='pass'):
                        qc_germ='pass'


                        if rt_t >75:
                            comb='hom'

                ##################################
        dic['qc'] = 'pass' if 'fail' not in qc_rec else 'fail'
        dic['somatic']='pass'
        if dic['qc']=='fail' and qc_germ=='pass':
            # print "2222222", dic,qc_germ
            dic['somatic']='fail'

            dic['qc']='pass'
            
            # print "test:",dic

            # print info
        frt = ':'.join([info['all_depth'], str(float(info['tumor_ratio%']) / 100), info['tumor_alt_forward_reverse']])
        #if dic['somatic']=='fail' and dic['qc']=='pass':## germline
        frt_normal = ':'.join([info['control_all_depth'],str(float(info['normal_ratio%']) / 100),info['normal_alt_forward_reverse']])
        # vcf_inf='DP='+info['all_depth']+';'+'AF='+str(float(info['tumor_ratio%'])/100)
        dic['comb']=comb
        dic['inf']=frt
        dic['norm_inf']=frt_normal
        #print qc_germ
        return dic
    else:
        ##add germline
        dp_cut = 20
        rt_somatic = 1
        dic = {}
        qc_rec = {}
        info = infos
        # print info
        # GT:DP:VD:ALD
        for i in info:
            # print i
            if onlymutect2: ##only mutect2 considered
                if i=='Software':
                    qc = 'pass' if info['Software'].find('mutect2') >=0 else 'fail'
                    qc_rec[qc]=1
            else:
                if i=='soft_count':###still use mutect2 as qc's software
                    if int(info[i])>=2:
                        qc='pass'
                        # print info['Software'].find('mutect2')
                    #     qc='pass' if info['Software'].find('mutect2') >=0 else 'fail'
                    else:
                        qc='fail'
                    qc_rec[qc]=1
                #print 1,qc
            if i=='all_depth' :
                # if info[i] == "." : info[i] = 0
                qc='pass' if int(info[i])>=dp_cut else 'fail'
                dp_germ='pass' if int(info[i])>=dp_cut else 'fail'
                qc_rec[qc] = 1
                #print 2,qc
            if i == 'tumor_ratio%':
                # if info[i] == ".": info[i] = 0
                qc = 'pass' if (float(info[i]) >= rt_somatic)  and \
                               (
                        int(re.split("[:,]",info["tumor_alt_forward_reverse"])[0])>=4 and
                        int(re.split("[:,]",info["tumor_alt_forward_reverse"])[1])>=2 and
                        int(re.split("[:,]",info["tumor_alt_forward_reverse"])[2])>=2
                ) else 'fail'
                qc_rec[qc] = 1
                # print 3,qc
            if i == "normal_ratio%":
                # if float(info[i]) == 0:
                    pass
                    # if i==
                # qc_rec[qc] = 1

                # if float(info[i])!=0 and float(info["tumor_ratio%"]) >=1 and float(float(info[i])/float(info["tumor_ratio%"])>=10)  and (int(re.split("[:,]",info["tumor_alt_forward_reverse"])[0])>=6and int(re.split("[:,]",info["tumor_alt_forward_reverse"])[1])>=3 and int(re.split("[:,]",info["tumor_alt_forward_reverse"])[2])>=3):
                #     qc="pass"
                # else:
                #     qc="fail"
                # qc_rec[qc] = 1
                # print qc_rec[qc]
                #print 4,qc
            if i=='tumor_alt_forward_reverse':
                rd=re.split('[:,]',info[i])
                # print rd
                qc='pass' if (
                        0.1<(float(rd[1])/float(rd[0]))<0.9
                        # int(rd[1])>0 and int(rd[2])>0
                              ) else 'fail'
                qc_rec[qc] = 1
                #print 5,qc
            if i in ['ExAC_Freq','1000G_ALL','ESP6500si_ALL','1000G_EAS','ExAC_EAS']:
                tag=poly_data(i,info[i])
                qc_rec[tag]=1
        # print qc_rec
        qc_germ = 'fail'
        comb = 'het'
        # if int(info["soft_count"]) >= 2:
        # # if info['Software'].find('vardictjava') >= 0:
        # #     if int(info['all_depth']) >= 20  :
        #         rt_t = float(info['tumor_ratio%'])
        #         # rt_s = float(info['normal_ratio%'])
        #         poly = 'pass'
        #         for i in ['ExAC_Freq', '1000G_ALL', 'ESP6500si_ALL', '1000G_EAS', 'ExAC_EAS']:
        #             tag = poly_data(i, info[i])
        #             # print i,info[i],tag
        #             if tag == 'fail':
        #                 poly = 'fail'
        #                 break
        #         cos = in_cosmic_times(info['COSMIC'])
        #         # print poly
        #         # print cos,rt_t,rt_s
        #         if rt_t >= 25 :
        #             if poly == 'pass' or (poly == 'fail' and cos == 'pass'):
        #                 qc_germ = 'pass'
        #                 if rt_t > 75: comb = 'hom'

        dic['qc'] = 'pass' if 'fail' not in qc_rec else 'fail'
        # print dic["qc"]
        dic['somatic'] = 'pass'
        if dic['qc'] == 'fail' and qc_germ == 'pass':
            dic['somatic'] = 'fail'
            dic['qc'] = 'pass'
        frt = ':'.join([info['all_depth'], str(float(info['tumor_ratio%']) / 100), info['tumor_alt_forward_reverse']])
        # if dic['somatic']=='fail' and dic['qc']=='pass':## germline
        frt_normal = ':'.join(
            [info['all_depth'], str(float(info['tumor_ratio%']) / 100), info['tumor_alt_forward_reverse']])
        # vcf_inf='DP='+info['all_depth']+';'+'AF='+str(float(info['tumor_ratio%'])/100)
        dic['comb'] = comb
        dic['inf'] = frt
        dic['norm_inf'] = frt_normal
        # print qc_germ
        # print dic
        return dic

def genetic_gene(f):
    ss=set()
    with open(f) as ff:
        for l in ff:
            ss.add(l.strip())
    return ss

def genetic_pos(f):
    ss=set()
    with open(f) as ff:
        for l in ff:
            field=l.strip().split('\t')
            ss.add('\t'.join(field[1:3]+['.']+field[3:5]))
    return ss


def main():
    parse = argparse.ArgumentParser('need vcf-sort software installed in the inviroment\n')
    parse.add_argument("--infile",help="input file")
    parse.add_argument("--code",help="sample code")
    parse.add_argument("--genetic",action='store_true',help="bool value, whether or not output genetic mutations")
    parse.add_argument("--genetic_filter",action='store_true',help="bool value, whether or not filter genetic mutations by genes and mutation poses")
    parse.add_argument("--genetic_gene",help="genetic gene list",default='/gpfs/software/pipeline/pipeline/database/local_database/gene_list/genetic_gene.txt')
    parse.add_argument("--genetic_pos",help="genetic pos list",default='/gpfs/software/pipeline/pipeline/database/local_database/genetic_db/518gene_3M.v3.txt')
    parse.add_argument("--sofstrict",action='store_true',help="whether or not restrict to mutect2, otherwise restrict to 3 software shared")
    parse.add_argument("--outfile",help="output file",default='for_watson.vcf')
    parse.add_argument("--pair",action='store_true', help= "pair yes/no")
    args=parse.parse_args()
    # print args.pair
    genetic_g=genetic_gene(args.genetic_gene)
    genetic_p=genetic_pos(args.genetic_pos)
    
    anno=glob.glob(args.infile)[0]
    try:
        mutect2=glob.glob(os.path.join(os.path.dirname( os.path.dirname(anno)),'snv/tumor/vardictjava/*.vardictjava.vcf'))[0]#'snv/pair/mutect2/*mutect2_paired.vcf'))[0]
        vcf_head_sample=os.popen('less '+mutect2+'|head -1000|awk \'/^##SAMPLE/\'').read()
    except IndexError:
        vcf_head_sample="##SAMPLE=<ID=NORMAL,File=-,SampleName=-\nSAMPLE=<ID=TUMOR,File=-,SampleName=-\n";

    if args.code:
        code=args.code
    else:
        code=os.popen('less '+anno+'|head -2|tail -1|awk -F"\t" \'{print $1}\'').read()
    #hd2='\t'.join(['#CHROM',  'POS',     'ID',      'REF',     'ALT',     'QUAL',    'FILTER',  'INFO',    'FORMAT',  code])
    hd2='\t'.join(['#CHROM',  'POS',     'ID',      'REF',     'ALT',     'QUAL',    'FILTER',  'INFO',    'FORMAT',  'NORMAL','TUMOR'])
    f = gzip.open(anno) if anno.endswith('.gz') else open(anno)

    o1=open(args.outfile,'w')
    o1tep=open(args.outfile+'.tep','w')
    o2=open(args.outfile+'.mut','w')
    vcf_head,vcf_head_contig=header()
    #print vcf_head,vcf_head_contig
    o1tep.write(vcf_head+'\n')
    o1tep.write(vcf_head_sample)
    #time.sleep( 5 )
    o1tep.write(vcf_head_contig+'\n')
    o1tep.write(hd2+'\n')

    head={}
    uniq={}
    for line in f:
        # print f
        info = {}
        fileds=line.split('\t')
        if line.startswith('Samplecode'):
            o2.write(line)
            for k,v in enumerate(fileds):
                if v=='MetaLR_pred':
                    break
                head[v]=k
                # print head[v]
            continue
        tg='\t'.join([fileds[head['chr']],fileds[head['pos']],'.',fileds[head['ref']],fileds[head['alt']]])
        for i in head:
            try:
                info[i] =fileds[head[i]]
            except IndexError:
                print "error",i,head[i]
        qc=QC_tissue(pair=args.pair,infos=info,onlymutect2=args.sofstrict)
        # print qc['qc']=='pass'
        if qc['qc']=='pass':
            ##if 3 software share while only mutect2 pass, this pos will be output still
            if tg in uniq:
                continue
            else:
                uniq[tg]=1
            comb = '0/1' if qc['comb'] == 'het' else '1/1'
            if qc['somatic']=='pass':
                tp = 'WfG_variant_origin=somatic'
                norm_comb='0/0'
            else:
                if args.genetic:
                    if args.genetic_filter :
                        if not( tg in genetic_p or info['GENE'] in genetic_g): continue
                
                    tp='WfG_variant_origin=germline'
                    norm_comb =comb
                else:
                    continue
            
            o2.write(line)

            tep='\t'.join([tg,'.','PASS',tp,'GT:DP:AF:VD:ALD',norm_comb+':'+qc['norm_inf'],comb+':'+qc['inf']])
            # print tep
            o1tep.write(tep+'\n')
            # print oltep
        # else:
        #     if info[i]


    o1tep.close()
    if os.path.exists(args.outfile):
        os.system('rm '+args.outfile)
        os.system('vcf-sort '+args.outfile+'.tep >'+args.outfile)
        os.system('rm '+args.outfile+'.tep');
    f.close()
    o2.close()
    

if __name__=='__main__':
    main()

