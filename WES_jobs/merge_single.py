#-*- coding:utf-8 -*-
import pandas as pd
import sys,os,glob,re
import argparse
import numpy as np
import subprocess
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument("-input",dest="input",help=u"输入需要分析文件")
parser.add_argument("-sample",dest="sample",help=u"样本编号 例子：/gpfs/home/zhaohongqiang/WES_linshi/merge.py -input PM2018091523.sofstric.watson.vcf -out PM2018091523.sofstric.watson.vcf -sample WES_T18001129")
args=parser.parse_args()
def makehash():
    return defaultdict(makehash)
Dicter=makehash()


def annover():##注释vcf脚本
    vcflist=glob.glob(os.path.abspath(args.input)+"/*.watson.vcf")[0]
    vcf_input=vcflist
    out_dir_file=vcflist
    annovar = "perl /gpfs/home/zhaohongqiang/pipeline/pipe_bin/table_annovar.pl {vcf_input} /gpfs/software/pipeline/pipeline/database/annovar_humandb -buildver hg19 -remove -protocol refGene,cosmic70,avsnp147,esp6500siv2_all,exac03,1000g2015aug_all,1000g2015aug_eas,clinvar_20160302,dbnsfp30a -operation g,f,f,f,f,f,f,f,f -nastring . -vcfinput -out {out_dir_file}".format(**locals())
    p=subprocess.Popen(annovar, shell=True)
    a=p.wait()
    print "进程状态：",a
    return out_dir_file



def Snpeff_run():
    vcflist = glob.glob(os.path.abspath(args.input) + "/*.watson.vcf")[0]
    vcf_input = vcflist
    out_dir_file = vcflist
    snpeff = "java -Xmx4G -jar /gpfs/home/zhaohongqiang/software/snpEff/snpEff.jar -c /gpfs/home/zhaohongqiang/software/snpEff/snpEff.config hg19 {vcf_input} > {out_dir_file}_snpeff.vcf".format(
        **locals())
    print snpeff
    p = subprocess.Popen(snpeff, shell=True)
    a = p.wait()
    print "进程状态：",a
    snpfile = "{out_dir_file}_snpeff.vcf".format(**locals())
    return snpfile



def snpeffdata(NM,cds,snpfile):

    a = pd.read_table(snpfile, sep="\t",skiprows=37,header=0)
    for i in a.INFO.str.split(";", expand=True)[1]:  # .str.split(",",expand=True)
        if [e for e in i.split(",") if e.find(NM) != -1 and e.find(cds) != -1]:
            sample = [e for e in i.split(",") if e.find(NM) != -1 and e.find(cds) != -1][0]
            if sample.split("|")[10]:

                AA = sample.split("|")[10]
            else:
                AA =np.nan
            print sample
            return AA
###'''处理注释后的# vcf'''###



files2=pd.read_table("/gpfs/home/zhaohongqiang/WES_linshi/WfG_canonical_transcript.txt")
out_dir_file=annover()
files_vcf=pd.read_table(out_dir_file+".hg19_multianno.vcf",sep="\t",skiprows=91,header=0)
files_INFO=files_vcf.INFO.str.split(";",expand=True)
files=pd.DataFrame(columns=[u'Chr', u'Start', u'End', u'Ref', u'Alt', u'Func.refGene',
       u'Gene.refGene', u'GeneDetail.refGene', u'ExonicFunc.refGene',
       u'AAChange.refGene', u'cosmic70', u'avsnp147', u'esp6500siv2_all',
       u'ExAC_ALL', u'ExAC_AFR', u'ExAC_AMR', u'ExAC_EAS', u'ExAC_FIN',
       u'ExAC_NFE', u'ExAC_OTH', u'ExAC_SAS', u'1000g2015aug_all',
       u'1000g2015aug_eas', u'CLINSIG', u'CLNDBN', u'CLNACC', u'CLNDSDB',
       u'CLNDSDBID', u'SIFT_score', u'SIFT_pred', u'Polyphen2_HDIV_score',
       u'Polyphen2_HDIV_pred', u'Polyphen2_HVAR_score', u'Polyphen2_HVAR_pred',
       u'LRT_score', u'LRT_pred', u'MutationTaster_score',
       u'MutationTaster_pred', u'MutationAssessor_score',
       u'MutationAssessor_pred', u'FATHMM_score', u'FATHMM_pred',
       u'PROVEAN_score', u'PROVEAN_pred', u'VEST3_score', u'CADD_raw',
       u'CADD_phred', u'DANN_score', u'fathmm-MKL_coding_score',
       u'fathmm-MKL_coding_pred', u'MetaSVM_score', u'MetaSVM_pred',
       u'MetaLR_score', u'MetaLR_pred', u'integrated_fitCons_score',
       u'integrated_confidence_value', u'GERP++_RS', u'phyloP7way_vertebrate',
       u'phyloP20way_mammalian', u'phastCons7way_vertebrate',
       u'phastCons20way_mammalian', u'SiPhy_29way_logOdds', u'Otherinfo',63,64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,76])

files["Func.refGene"]=files_INFO[2].str.split("=",expand=True)[1]
files["ExonicFunc.refGene"]=files_INFO[5].str.split("=",expand=True)[1]
files[72]=files_INFO[0]
files["Gene.refGene"]=files_INFO[3].str.split("=",expand=True)[1]
files[75]=files_vcf.TUMOR
files[74]=files_vcf.NORMAL
files["Ref"]=files_vcf.REF
files["Alt"]=files_vcf.ALT
files["Chr"]=files_vcf["#CHROM"]
files["Start"]=files_vcf.POS
files["GeneDetail.refGene"]=files_INFO[4].str.split("=",expand=True)[1]
files["AAChange.refGene"]=files_INFO[6].str.split("=",expand=True)[1]
files["cosmic70"]=files_INFO[7].str.split("=",expand=True)[1].astype(str).str.replace(r"\Sx3d","=").str.replace(r"\Sx3b",";")
################################################################################################################################################
# germlinelist=files[(files["Func.refGene"]!="downstream")&(files["Func.refGene"]!="intergenic")&(files["Func.refGene"]!="intronic")&(files["Func.refGene"]!="ncRNA_intronic")&(files["Func.refGene"]!="upstream")&(files["Func.refGene"]!="UTR3")&(files["Func.refGene"]!="UTR5")&(files["ExonicFunc.refGene"]!="synonymous SNV")&("germline" == files[72].str.split("=",expand=True)[1])]
somaticlist=files[
    (files["Func.refGene"]!="downstream")&
    (files["Func.refGene"]!="intergenic")&
    (files["Func.refGene"]!="intronic")&
    (files["Func.refGene"]!="ncRNA_intronic")&
    (files["Func.refGene"]!="upstream")&
    (files["Func.refGene"]!="UTR3")&
    (files["Func.refGene"]!="UTR5")&
    (files["ExonicFunc.refGene"]!="synonymous_SNV")
]
##############################################################################################################################################################################################################################################################################################################################################
annofile=glob.glob("{}/*.txt_ti_all_leve.xls".format(os.path.abspath(args.input)))[0]#'/gpfs/project/RD_bioinfor/WES_clinical/20180926_WES/wangtong_ffpe/annotation/reformed.PM2018091523.3tool.hg19_multianno.txt_ti_all_leve.xls'
mutfile = glob.glob("{}/*.vcf.mut".format(os.path.abspath(args.input)))[0]
mut=pd.read_table(mutfile,sep="\t",header=0,encoding = 'gbk')
level=pd.read_table(annofile,sep="\t",header=0,encoding = 'gbk')




def to_merge(files,ctrl=1):
    df=pd.DataFrame()
    columns = ["Samplecode", "tag", 'Class', 'Reported', 'Score', 'GENE532', 'Software', 'GENE', 'CDS_Mutation',
                  'AA_Mutation', 'Exon', 'tumor_ratio', 'normal_ratio', 'Leve', 'QC', 'Likely_germline', 'NM', 'status',
                  'target_chemo_mark', 'COSMIC', 'COSMIC_somatic_status', 'ExonicFunc', 'alt_depth', 'all_depth',
                  'control_alt_depth', 'control_all_depth', 'tumor_alt_forward_reverse', 'normal_alt_forward_reverse',
                  'homopolymer', 'dbsnp', 'ALL_1000G', 'EAS_1000G', 'ExAC_Freq', 'ExAC_EAS', 'ESP6500si_ALL',
                  'ClinVar_SIG', 'ClinVar_DIS', 'Pathogenicity', 'soft_count', 'Freebayes_num', 'Vardict_num',
                  'Mutation_db', 'DISEASE', 'FDACFDA1', 'FDACFDA2', 'DRUGS', 'INFORMATION', 'COSMIC_ID', 'chr', 'pos',
                  'gene_symbol', 'ref', 'alt', 'SIFT_pred', 'MutationAssessor_pred', 'Polyphen2_HDIV_pred',
                  'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'FATHMM_pred', 'PROVEAN_pred',
                  'fathmm_MKL_coding_pred', 'MetaSVM_pred', 'MetaLR_pred', 'info1', 'info2', 'info3', 'info4', 'info5',
                  'info6', 'info7', 'info8', 'info9', 'info10', 'Chr', 'Position', 'Ref', 'A', 'T', 'C', 'G', 'N']
    df["gene"]=files["Gene.refGene"]
    df["ExonicFunc"]=files["ExonicFunc.refGene"]
    df["tumor_ratio%"]=files[75][files[75]!="."].str.split(":",expand=True)[2].astype(float)*100
    df["nomal_ratio%"]=files[74][files[74]!="."].str.split(":",expand=True)[2].astype(float)*100
    df["ref"]=files["Ref"]
    df['alt']=files["Alt"]
    df["pos"]=files["Start"]
    df["chr"]=files["Chr"]
    df["COSMIC"]=files["cosmic70"]
    df["COSMIC_ID"] = files["cosmic70"][files["cosmic70"] != "."].str.split("[;,]", expand=True)[0]
    df=df.reset_index(drop=True)



    #############################匹配wason##################################
    def get_tp_anno(tp_anno):
                    all_tp=[x.split(":") for x in tp_anno.split(",")] #[[gene,tp,exon,c_anno,(p_anno)],[]]
                    if len(all_tp[0])>=4 :

                            # print all_tp
                            c_anno=[x[3] for x in all_tp if len(x)>=4]
                            tp=[x[1] for x in all_tp]
                            genename=[x[0] for x in all_tp]
                            exon=[x[2] for x in all_tp]
                            p_anno=[x[4] for x in all_tp if len(x)>=5]
                    else:
                            c_anno,p_anno,tp,exon=[],[],[],[]
                    return c_anno,p_anno,tp,exon

    def get_splicing(tp_anno):
                    all_tp = [x.split(":") for x in tp_anno.split(",")]  # [[gene,tp,exon,c_anno,(p_anno)],[]]
                    if len(all_tp[0]) >= 3:

                        print all_tp
                        c_anno = [x[2] for x in all_tp]
                        tp = [x[0] for x in all_tp]
                        # genename = [x[0] for x in all_tp]
                        exon = [x[1] for x in all_tp]
                        p_anno =[] #[x[3] for x in all_tp if len(x) >= 5]
                    else:
                        c_anno, p_anno, tp, exon = [], [], [], []
                    return c_anno,p_anno, tp, exon
    NM=[]
    nums=0
    sm=0
    bums=0
    genelist=[]
    plist=[]
    clist=[]
    exonlist=[]
    for g,i,gene in zip(files["GeneDetail.refGene"],files["AAChange.refGene"],files["Gene.refGene"]):
        sm+=1
        if i!=".":
            c, p, tp,exon= get_tp_anno(i)
        else:
            # print g
            c, p, tp, exon = get_splicing(g)
        # print tp
        if tp == []:
            bums += 1
            exonlist.append(np.nan)
            NM.append(np.nan)
            genelist.append(np.nan)
            plist.append(np.nan)
            clist.append(np.nan)
        else:
            exonlist.append(np.nan)
            genelist.append(np.nan)
            clist.append(np.nan)
            plist.append(np.nan)
            NM.append(np.nan)
            for geneses, nmber in zip(files2['Gene_symbol'].str.split(".", expand=True)[0],
                                        files2['primary_transcript.refseqID'].str.split(".", expand=True)[0]):
                    if nmber in tp :
                        nums += 1
                        if i != ".":
                            AA = pd.DataFrame(i.split(","))
                            AA = AA[AA[0].str.split(":", expand=True)[1] == nmber]
                            SS = AA[0].str.split(":", expand=True)
                            # print SS
                            if 4 in SS.columns:
                                listp = [pk for pk in AA[0].str.split(":", expand=True)[4]][0]
                            else:
                                listp = np.nan
                            if 3 in SS.columns:
                                listc = [ck for ck in AA[0].str.split(":", expand=True)[3]][0]
                            else:
                                listc = np.nan
                            listexom = [ek for ek in AA[0].str.split(":", expand=True)[2]][0]
                        else:

                            AA = pd.DataFrame(g.split(","))
                            print AA
                            AA = AA[AA[0].str.split(":", expand=True)[0] == nmber]
                            # print AA
                            SS = AA[0].str.split(":", expand=True)

                            listp = np.nan
                            if 2 in SS.columns:
                                listc = [ck for ck in AA[0].str.split(":", expand=True)[2]][0]
                            else:
                                listc = np.nan
                            listexom = [ek for ek in AA[0].str.split(":", expand=True)[1]][0]

                        # print "wason 转录本：",nmber,"vcf_annover:",tp,p,"out",listp,listc
                        del genelist[sm-1]
                        del plist[sm-1]
                        del clist[sm-1]
                        del NM[sm-1]
                        del exonlist[sm-1]
                        genelist.insert(sm,geneses)
                        exonlist.insert(sm, listexom)
                        plist.insert(sm, listp)
                        clist.insert(sm,listc)
                        NM.insert(sm,nmber)


    # print len(plist)
    # print plist
    df["genelist"]=pd.DataFrame(genelist)
    df["AA_mutation"]=pd.DataFrame(plist)
    df["exon"]=pd.DataFrame(exonlist)
    df["CDS_Mutation"]=pd.DataFrame(clist)
    df["NM"]=pd.DataFrame(NM)
    ###############################突变匹配wason###########################################
    df = df[df["NM"] == df["NM"]]
    df = df.drop_duplicates(["CDS_Mutation","AA_mutation","NM"])
    df=pd.merge(df, mut[["CDS_Mutation", "AA_Mutation", 'GENE', "chr", "pos", "ref", "alt"]],
             on=["chr", "pos", "ref", "alt"], how="left")
    df["AA_mutation"] = df["AA_mutation"].fillna(df[df.AA_mutation!=df.AA_mutation]["AA_Mutation"]).replace("p.(=)",np.NaN)
    df["CDS_Mutation"] = df["CDS_Mutation_x"]
    def gzdata():
        sc = df["CDS_Mutation"].str.extract(r"c.([AGCT]+)([0-9]+)([AGCT]+)")
        return "c." + sc[1] + sc[0] + ">" + sc[2]

    df["CDS_Mutation"]=gzdata().fillna(df["CDS_Mutation"][df["CDS_Mutation"].str.extract(r"c.([AGCT]+[0-9]+[AGCT]+)")[0].isna()])
    # df["CDS_Mutation"]=df["CDS_Mutation_x"]
    df.drop("CDS_Mutation_y",axis=1)
    df["AA_mutation"]=df["AA_mutation"].replace(".",np.nan)
    # print df["AA_mutation"][df["AA_mutation"]!=df["AA_mutation"]]
    snpfile = Snpeff_run()
    for aa, NM,pos,cds  in zip(
            df["AA_mutation"][df["AA_mutation"]!=df["AA_mutation"]].index,
            df["NM"][df["AA_mutation"]!=df["AA_mutation"]],
            df["pos"][df["AA_mutation"]!=df["AA_mutation"]],
            df["CDS_Mutation"][df["AA_mutation"]!=df["AA_mutation"]]
    ):
        cds = re.findall(r"(?=.)\d+", cds)[0]
        print NM, aa,cds
        AA=snpeffdata(NM,cds,snpfile)
        df.loc[aa,"AA_mutation"]=AA
    print df
    df.reset_index(drop=True)
    return df


def Somatic_positive():
        import pymysql
        conn = pymysql.connect(host='10.2.22.24', port=3306, user='bioinford', passwd='223355', db='vcangene')
        sql = "select * from seq_benign_somatic"
        sp = pd.read_sql(sql, conn)

        gmlinejd["ref"] = gmlinejd["Ref"]
        gmlinejd["alt"] = gmlinejd["Alt"]
        gmlinejd["pos"] = gmlinejd["Start"]
        for chr, ref, alt, pos in zip(gmlinejd.chr, gmlinejd.ref, gmlinejd.alt, gmlinejd.pos):
            # print chr , ref, alt ,pos
            if [(pt.chr == chr) & (pt.ref == ref) & (pt.alt == alt) & (pt.pos == pos)] == True:
                pt = pt.drop(pt.index)


def Gearmline_positive():
        import pymysql
        conn = pymysql.connect(host='10.2.22.24', port=3306, user='bioinford', passwd='223355', db='vcangene')
        sql = "select * from seq_benign_germline"
        gp = pd.read_sql(sql, conn)
        return gp


if __name__=="__main__":
    # df=to_merge(germlinelist,ctrl=1)
    dr=to_merge(somaticlist,ctrl=2)
    ds=pd.DataFrame(columns=["Samplecode","tag",'Class','Reported','Score', 'GENE532', 'Software', 'GENE', 'CDS_Mutation', 'AA_Mutation', 'Exon', 'tumor_ratio', 'normal_ratio', 'Leve', 'QC', 'Likely_germline', 'NM', 'status', 'target_chemo_mark', 'COSMIC', 'COSMIC_somatic_status', 'ExonicFunc', 'alt_depth', 'all_depth', 'control_alt_depth', 'control_all_depth', 'tumor_alt_forward_reverse', 'normal_alt_forward_reverse', 'homopolymer', 'dbsnp', 'ALL_1000G', 'EAS_1000G', 'ExAC_Freq', 'ExAC_EAS', 'ESP6500si_ALL','ClinVar_SIG', 'ClinVar_DIS','Pathogenicity','soft_count', 'Freebayes_num', 'Vardict_num','Mutation_db', 'DISEASE', 'FDACFDA1', 'FDACFDA2', 'DRUGS', 'INFORMATION','COSMIC_ID', 'chr', 'pos', 'gene_symbol', 'ref', 'alt', 'SIFT_pred', 'MutationAssessor_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'FATHMM_pred', 'PROVEAN_pred', 'fathmm_MKL_coding_pred', 'MetaSVM_pred', 'MetaLR_pred'])
    dg=pd.DataFrame(columns=['Samplecode', 'taged', 'Class', 'Reported', 'Score', '532GENE',
           'Software', 'GENE', 'CDS_Mutation', 'AA_Mutation', 'Exon',
           'normal_ratio%', 'tumor_ratio%', 'Leve', 'QC', 'Likely_germline',
           'NM', 'status', 'target_chemo_mark', 'COSMIC',
           'COSMIC_somatic_status', 'ExonicFunc', 'alt_depth', 'all_depth',
           'control_alt_depth', 'control_all_depth',
           'tumor_alt_forward_reverse', 'normal_alt_forward_reverse',
           'homopolymer', 'dbsnp', '1000G_ALL', '1000G_EAS', 'ExAC_Freq',
           'ExAC_EAS', 'ESP6500si_ALL', 'ClinVar_SIG', 'ClinVar_DIS',
           'Pathogenicity', 'soft_count', 'Freebayes_num', 'Vardict_num',
           'Mutation_db', 'DISEASE', 'FDACFDA1', 'FDACFDA2', 'DRUGS',
           'INFORMATION', 'COSMIC_ID', 'chr', 'pos', 'gene_symbol', 'ref',
           'alt', 'SIFT_pred', 'MutationAssessor_pred', 'Polyphen2_HDIV_pred',
           'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred',
           'FATHMM_pred', 'PROVEAN_pred', 'fathmm-MKL_coding_pred',
           'MetaSVM_pred', 'MetaLR_pred' ])
    PM=args.sample
    ds['GENE']=dr["genelist"]
    ds["CDS_Mutation"]=dr["CDS_Mutation"]
    ds["AA_Mutation"]=dr["AA_mutation"]
    ds["Exon"]=dr["exon"]
    ds["NM"]=dr["NM"]
    ds["ref"]=dr["ref"]
    ds["alt"]=dr["alt"]
    ds["pos"]=dr["pos"]
    ds["chr"]=dr["chr"]
    ds['COSMIC']=dr["COSMIC"]
    ds['tumor_ratio']=dr['tumor_ratio%']
    ds['normal_ratio']=dr['nomal_ratio%']
    ds["COSMIC_ID"]=dr["COSMIC_ID"]
    ds["Samplecode"]=ds["Samplecode"].fillna("{}".format(PM))

    # dg['GENE']=df["genelist"]
    # dg["CDS_Mutation"]=df["CDS_Mutation"]
    # dg["AA_Mutation"]=df["AA_mutation"]
    # dg["Exon"]=df["exon"]
    # dg["NM"]=df["NM"]
    # dg["ref"]=df["ref"]
    # dg["alt"]=df["alt"]
    # dg["pos"]=df["pos"]
    # dg["chr"]=df["chr"]
    # dg['tumor_ratio%']=df['tumor_ratio%']
    # dg['normal_ratio%']=df['nomal_ratio%']
    # dg['COSMIC']=df["COSMIC"]
    # dg["COSMIC_ID"]=df["COSMIC_ID"]
    # dg["Samplecode"]=dg["Samplecode"].fillna("{}".format(PM))

    ds.to_excel("{}/ALL_{}_somatic.xls".format(args.input,PM),index=False)
    # dg.to_excel("{}/ALL_{}_genetic.xls".format(args.input,PM),index=False)


    ###########合并allevel##############
    # pt=pd.merge(dg[['GENE',"chr","pos","ref","alt","NM","CDS_Mutation","AA_Mutation"]],level,on=["chr","pos","ref","alt"])
    # pt=pt[pt.Software=="vardictjava-pair"].drop_duplicates(["chr","pos","ref","alt"])
    pso=pd.merge(ds[['GENE',"chr","pos","ref","alt","NM","CDS_Mutation","AA_Mutation"]],mut,on=["chr","pos","ref","alt",'NM'])
    b=pso[[
        u'chr',
        u'pos', u'ref',
        u'alt', u'NM',
        u'Software', u'GENE_x'
    ]]["Software"]
    pso=pso[(b.isin(["vardictjava-single"])) | (b.isin(["freebayes-single"]))].drop_duplicates(
        ["chr", "pos", "ref", "alt", 'NM'], keep="last")
    # pso=pso[pso.Software=="mutect2-pair"].drop_duplicates(["chr","pos","ref","alt"])

    #########合并###############






    colums=[u'Samplecode', u'taged', u'Class', u'Reported', u'Score', u'532GENE',
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
    colums2 = [                       u'chr',                        u'pos',
                                  u'ref',                        u'alt',

                                         u'Samplecode',
                                u'taged',                      u'Class',
                             u'Reported',                      u'Score',
                              u'532GENE',                   u'Software',
                                 u'GENE',
                                              u'Exon',
                         u'tumor_ratio%',              u'normal_ratio%',
                                 u'Leve',                         u'QC',
                      u'Likely_germline',
                               u'status',          u'target_chemo_mark',
                               u'COSMIC',      u'COSMIC_somatic_status',
                           u'ExonicFunc',                  u'alt_depth',
                            u'all_depth',          u'control_alt_depth',
                    u'control_all_depth',  u'tumor_alt_forward_reverse',
           u'normal_alt_forward_reverse',                u'homopolymer',
                                u'dbsnp',                  u'1000G_ALL',
                            u'1000G_EAS',                  u'ExAC_Freq',
                             u'ExAC_EAS',              u'ESP6500si_ALL',
                          u'ClinVar_SIG',                u'ClinVar_DIS',
                        u'Pathogenicity',                 u'soft_count',
                        u'Freebayes_num',                u'Vardict_num',
                          u'Mutation_db',                    u'DISEASE',
                             u'FDACFDA1',                   u'FDACFDA2',
                                u'DRUGS',                u'INFORMATION',
                            u'COSMIC_ID',                u'gene_symbol',
                            u'SIFT_pred',      u'MutationAssessor_pred',
                  u'Polyphen2_HDIV_pred',        u'Polyphen2_HVAR_pred',
                             u'LRT_pred',        u'MutationTaster_pred',
                          u'FATHMM_pred',               u'PROVEAN_pred',
               u'fathmm-MKL_coding_pred',               u'MetaSVM_pred',
                          u'MetaLR_pred',                       u'是否报出',
                                 u'是否有药',                         u'功能',
                                 u'批准疾病',                u'用于批准疾病的敏感用药',
                          u'用于其他疾病的敏感用药',                     u'临床试验用药',
                               u'其他耐药信息',                         u'级别',
                                   u'注释',                        u'Chr',
                             u'Position',                        u'Ref',
                                    u'A',                          u'T',
                                    u'C',                          u'G',
                                    u'N']
    # pd.merge(dg[["chr","pos","ref","alt"]],level,on=["chr","pos","ref","alt"],how="left").drop_duplicates(["chr","pos","ref","alt"])
    # pt[colums2]=pt[colums2].replace(".",np.nan)
    # pt["othergene"]=pt["gene_symbol"]
    # pt["otherCDS"]=pt["CDS_Mutation_y"]
    # pt["otherAA"]=pt["AA_Mutation_y"]
    # pt["NM"]=pt["NM_x"]
    # pt['GENE']= pt['GENE_x']
    # pt["CDS_Mutation"]=pt["CDS_Mutation_x"]
    # pt["AA_Mutation"]=pt["AA_Mutation_x"]
    # pt=pt.drop(['GENE_x','GENE_y',"NM_x","NM_y","CDS_Mutation_x","CDS_Mutation_y","AA_Mutation_x","AA_Mutation_y"],axis=1)
    # pt[colums2] = pt[colums2].replace(".", np.nan)
    # ptlist=pt[colums]




    # pso["othergene"]=pso["gene_symbol"]
    # pso["otherCDS"]=pso["CDS_Mutation_y"]
    # pso["otherAA"]=pso["AA_Mutation_y"]
    pso["NM"]=pso["NM"]
    pso['GENE']=pso['GENE_x']
    pso["CDS_Mutation"]=pso["CDS_Mutation_x"]
    pso["AA_Mutation"]=pso["AA_Mutation_x"]
    pso=pso.drop(["CDS_Mutation_x","CDS_Mutation_y","AA_Mutation_x","AA_Mutation_y"],axis=1)
    pso[colums2] = pso[colums2].replace(".", np.nan)
    psolist=pso[colums]
    # ptlist.to_excel("{}/READ_{}_genetic.xls".format(args.input,PM),index=False,encoding="uft-8")
    psolist.to_excel("{}/READ_{}_somatic.xls".format(args.input,PM),index=False)



