#coding:utf-8
__author__ = "JC"
__date__ = '2018/11/30 0030 上午 9:40'
import MySQLdb
import os,json
import pandas as pd

drug_info_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'database', 'cancer_info.txt')

class Cancer_Info_From_lims:
    def __init__(self):
        self.db = MySQLdb.connect("10.2.22.26", "report", "bioreport123$", "vlisdb", unix_socket="/tmp/mysql.socket", charset='utf8')
        self.cursor = self.db.cursor()
        self.log = []
        self.cancer_info = self._get()
    def _get(self):
        sql = "SELECT u.fullName '销售', sl.programName '检测项目', ss.newCode '样本编号', c.`name` '受检人姓名', c.gender '性别', c.age '年龄'," \
              "c.identity '身份证号', o.submitDate '送检时间', o.source '来源', h.`name` '医院', dp.`name` '科室', d.`name` '主治医生', ss.`code` '条码号'," \
              "ss.`level` '样本等级', ss.specimenType '样本类型', ss.specimenTypeRemark '样本类型说明', ss.sampleSize '采样数量', o.collectDate '采样日期'," \
              "s.arriveDate '收样日期', ss.sampleUnit '样本规格', ss.expressType '快递类型', ss.courierNumber '快递单号', co.`name` '报告接收人',co.tel '报告接收电话'," \
              "co.address '报告接收地址', di.`name` '疾病种类', dc.`name` '分型', o.clinicalDiagnosis '临床诊断', o.primarySite '原发部位', o.transferSite '转移部位'," \
              "cl.`name` '临床分期', o.stagingDetail '分期详情', o.familyHistory '有无家族史', o.smokingHistory '有无吸烟史', o.quitSmoking '停止吸烟时长'," \
              "o.recentTest '最近检测基因', o.recentTestDate '最近检测时间', o.recentTestMethod '最近检测方法', o.recentTestResult '最近检测结果', o.remark1 '备注1' " \
              "FROM salessample ss LEFT JOIN sample s ON s.id = ss.sample_id LEFT JOIN samplelink sl ON ss.sampleLink_id = sl.id LEFT JOIN orders o ON ss.orders_id = o.id " \
              "LEFT JOIN hospital h ON o.hospital_id = h.id LEFT JOIN department dp ON o.department_id = dp.id LEFT JOIN `user` u ON o.sales_id = u.id " \
              "LEFT JOIN client c ON o.client_id = c.id LEFT JOIN doctor d ON o.attendingDoctor_id = d.id LEFT JOIN contact co ON o.reception_id = co.id " \
              "LEFT JOIN disease di ON o.disease_id = di.id LEFT JOIN diseaseclassification dc ON o.diseaseClass_id = dc.id LEFT JOIN clinicalstage cl ON o.clinicalStage_id = cl.id " \
              "WHERE ss.sampleLink_id IS NOT NULL;"
        self.cursor.execute(sql)
        result = self.cursor.fetchall()
        cancer_info = {}
        sample_code_dup = []
        for all in result:
            sample_code = all[2]
            cancer_category = all[25]
            if sample_code in sample_code_dup:
                self.log.append("%s 有可能有重复")
            sample_code_dup.append(sample_code)
            cancer_info[sample_code] = cancer_category

        return cancer_info

class Drug_Info:
    def __init__(self , drug_file):
        self.drug_file = drug_file
    def parse_file(self):
        import codecs
        result = {}
        with codecs.open(self.drug_file, encoding='gbk') as f:
            head = next(f).split("\t")
            result['heads'] = head
            for line in f:

                tmp_list = line.split("\t")

                gene = tmp_list[ head.index("GENE") ].strip()
                cds = tmp_list[ head.index("CDS_Mutation") ].strip()
                AA = tmp_list[ head.index("AA_Mutation") ].strip()
                exon = tmp_list[ head.index("Exon") ].strip()
                cancer = tmp_list[head.index("Score")]
                level = tmp_list[head.index(u"级别")]
                if exon == "exon19":
                    key = "_".join([cancer ,gene , cds  , exon])
                else:
                    key = "_".join([cancer , gene , cds , AA , exon])
                if key not in result:
                    result[key]= {}
                result[key][level] = tmp_list
        return result

class Mutation_Info:
    def __init__(self , mutation_file , drug_dict, cancer_info , out_file):
        self.mutation_file = mutation_file
        self.drug_dict = drug_dict
        self.cancer_info = cancer_info
        self.out_file =out_file
        if not os.path.isfile(self.mutation_file):
            raise Exception("%s not exists"%self.mutation_file)

    def _add_colorectal_cancer_wild_type(self, gene , headers , sample_code):
        a_list = ['']*len(headers)
        a_list[headers.index('Samplecode')] = sample_code
        a_list[headers.index('Score')] = u'结直肠癌'
        a_list[headers.index("GENE")] = gene
        a_list[headers.index("CDS_Mutation")] = u"野生型"
        return a_list

    def parse_file(self):
        df = pd.read_excel(self.mutation_file , dtype={"Samplecode": str} )
        df = df.fillna("")
        sample_code = df.ix[0, 'Samplecode']
        if sample_code in self.cancer_info:
            self.cancer = self.cancer_info[sample_code]
        else:
            self.cancer = "无"

        if self.cancer == u'结直肠癌':
            df.loc[df.shape[0]] = self._add_colorectal_cancer_wild_type('KRAS', list(df.columns) , sample_code)
            df.loc[df.shape[0]] = self._add_colorectal_cancer_wild_type('NRAS', list(df.columns) , sample_code)
            df.loc[df.shape[0]] = self._add_colorectal_cancer_wild_type('BRAF', list(df.columns), sample_code)


        result = {}
        result[0] = df.columns
        num =1
        # print df.columns
        for i in df.index:
            sample_code = df.ix[i, 'Samplecode']
            ###add cancer info from lims
            if sample_code in self.cancer_info:
                df.ix[i , 'Score'] = self.cancer_info[sample_code]
            else:
                df.ix[i, 'Score'] = u"没有从Lims中抓取到癌症信息"

            gene = df.ix[i , 'GENE']
            cds = df.ix[i, 'CDS_Mutation']
            AA = df.ix[i, 'AA_Mutation']
            exon = df.ix[i, 'Exon']
            cancer = df.ix[i, 'Score']
            ####CNV FUSION
            taged = df.ix[i , 'taged']
            if taged == u"CNV无" or taged == u'fusion无':
                result[num] = df.ix[i,]
                num+=1
                continue


            if exon == "exon19":
                key = "_".join([cancer ,gene, cds, exon])
            else:
                key = "_".join([cancer , gene, cds, AA, exon])


            if key in self.drug_dict:
                for level in self.drug_dict[key]:
                    drug_list = self.drug_dict[key][level]
                    drug_heads = self.drug_dict['heads']
                    df.ix[i, u"用于批准疾病的敏感用药"] = drug_list[drug_heads.index(u"用于批准疾病的敏感用药")].strip('"')
                    df.ix[i, u"用于其他疾病的敏感用药"] = drug_list[drug_heads.index(u"用于其他疾病的敏感用药")].strip('"')
                    df.ix[i, u"临床试验用药"] = drug_list[drug_heads.index(u"临床试验用药")].strip('"')
                    df.ix[i, u"其他耐药信息"] = drug_list[drug_heads.index(u"其他耐药信息")].strip('"')
                    df.ix[i, u"级别"] = drug_list[drug_heads.index(u"级别")].strip('"')
                    df.ix[i, u"注释"] = drug_list[drug_heads.index(u"注释")].strip('"')
                    result[num] = df.ix[i,]
                    num+=1
            else:
                result[num] = df.ix[i,]
                num +=1



        df_output = pd.DataFrame(result).T
        df_output.to_excel(self.out_file , index =False, header = False)


def main(infile , cancer_info, drug_dict , out_file):
    drug_info_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'database', 'cancer_info.txt')
    mutation_file = infile


    # drug_obj = Drug_Info(drug_info_file)
    # drug_dict = drug_obj.parse_file()

    # cancer_obj = Cancer_Info_From_lims()
    # cancer_info = cancer_obj.cancer_info

    mutation_obj = Mutation_Info(mutation_file, drug_dict,cancer_info ,out_file)
    mutation_obj.parse_file()


if __name__ == '__main__':
    drug_info_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'database', 'cancer_info.txt')
    mutation_file = ur"D:\zyxh\商业样本\13基因血液\突变位点\20181128\hufengjuan_ffpe\hufengjuan_ffpe_test.xls"
    out_file = os.path.join(os.path.dirname(mutation_file), 'test_test.xls')

    drug_obj = Drug_Info(drug_info_file)
    drug_dict = drug_obj.parse_file()

    cancer_obj = Cancer_Info_From_lims()
    cancer_info = cancer_obj.cancer_info

    mutation_obj = Mutation_Info(mutation_file, drug_dict,cancer_info ,out_file)
    mutation_obj.parse_file()

    # print json.dumps(drug_obj.parse_file(), indent=4)