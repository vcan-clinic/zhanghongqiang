# -*- coding: utf-8 -*-  
import os
from reportlab.pdfgen.canvas import Canvas  
from reportlab.pdfbase import pdfmetrics  
from reportlab.pdfbase.ttfonts import TTFont 
pdfmetrics.registerFont(TTFont('wqy', os.path.join(os.path.dirname(__file__),'wqy-microhei.ttc'))) 
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer,Image,Table,TableStyle
import time
import json
import sys

def rpt(json_report, pdf_report, sample_code, library_code, curr_date = ""):

    with open(json_report, 'r') as json_f:
        results_dict = json.load(json_f) 
    print(results_dict)

    def table_style(mut_list,label_column):
        trans = lambda x:x[:label_column]+[{"yes":u"检出","no":u"未检出"}[x[label_column]]]
        label_colour = []
        new_mut_list = []
        for n,m in enumerate(mut_list):
            new_mut_list.append(trans(m))
            if m[label_column] == "yes":
                label_colour.append(('BACKGROUND',(label_column,n+1),(label_column,n+1), colors.lightsalmon))
        return new_mut_list,label_colour
    results_dict["snvList"],snv_label_colour = table_style(results_dict["snvList"],-1)
    results_dict["fusionList"],fusion_label_colour = table_style(results_dict["fusionList"],-1)
 
    story=[]
    stylesheet=getSampleStyleSheet()
    normalStyle = stylesheet['Normal']
    if not curr_date:
       # curr_date = time.strftime("%Y年%m月%d日 %H:%M", time.localtime())
        curr_date = time.strftime("%Y年%m月%d日", time.localtime())   
    #标题：段落的用法详见reportlab-userguide.pdf中chapter 6 Paragraph
    rpt_title1 = '<para autoLeading="off" align="center"><b><font face="wqy" fontSize=14>人EGFR/ALK/BRAF/KRAS基因突变联合检测试剂盒(可逆末端终止测序法)</font></b><br/><br/></para>'
    story.append(Paragraph(rpt_title1,normalStyle)) 
    rpt_title2 = '<para autoLeading="off" align="center"><b><font face="wqy" fontSize=14>检测报告</font></b><br/><br/><br/></para>'
    story.append(Paragraph(rpt_title2,normalStyle))    

    sample_codes = '''<para autoLeading="off" align=left><b><font face="wqy" fontSize=12>样本编号</font></b><br/>
    <font face="wqy" fontSize=8>样本编号：%s</font><br/>
    <font face="wqy" fontSize=8>文库编号：%s</font><br/>
    <font face="wqy" fontSize=8>检测日期：%s</font><br/><br/></para>
    ''' %(sample_code, library_code, curr_date)
    story.append(Paragraph(sample_codes,normalStyle))

    type_code = "FFPE"
    if results_dict['qcFailed'] > 0 :
        control_code = "不满足要求"
    else:
        control_code = "满足要求" 
    sample_info = '''<para autoLeading="off" align=left><b><font face="wqy" fontSize=12>样本信息</font></b><br/>
    <font face="wqy" fontSize=8>样本类型：%s；样本质控情况：%s</font><br/><br/></para>
    ''' %(type_code, control_code)
    story.append(Paragraph(sample_info,normalStyle))
   
    if results_dict['snvDetected'] > 0 or results_dict['fusionDetected'] > 0:
        result_code = "检测结果中出现基因突变项"
    else:
        result_code = "检测结果中未出现基因突变项"
    test_info = '''<para autoLeading="off" align=left><b><font face="wqy" fontSize=12>检测提示</font></b><br/>
    <font face="wqy" fontSize=8>%s</font><br/><br/></para>
    ''' %(result_code)
    story.append(Paragraph(test_info,normalStyle))

    text = '''<para autoLeading="off" align=left><b><font face="wqy" fontSize=12>检测声明</font></b><br/>
    <font face="wqy" fontSize=8>1. 报告结果只包含 Panel 覆盖的基因位点的突变情况，突变形式为 DNA 的 SNV 和 InDel 突变及 RNA 的融合类型;</font><br/>
    <font face="wqy" fontSize=8>2. 报告结果只对本次受检样本负责，结果仅用于辅助癌症用药检测诊断，对于本次检测结果的解释，请进一步咨询医生;</font><br/>
    <font face="wqy" fontSize=8>3. 受检者需提供完整、准确、详实的个人资料。因受检者提供的用户资料不实造成其他误导因素而导致检测服务的中断、结果不准确，医院对此不承担责任。</font><br/><br/>
    </para>'''
    story.append(Paragraph(text,normalStyle))

    text = '<para autoLeading="off"><b><font face="wqy" fontSize=12>检测结果</font></b><br/><br/><font face="wqy" fontSize=8>1.点突变、小片段插入或缺失</font><br/><br/></para>'
    story.append(Paragraph(text,normalStyle))

    #图片，用法详见reportlab-userguide.pdf中chapter 9.3 Image

    #表格数据：用法详见reportlab-userguide.pdf中chapter 7 Table
    snv_head = ['基因', 'Cosmic ID', '碱基突变', '氨基酸突变', '是否检出'] 
    results_dict['snvList'].insert(0, snv_head)
    #创建表格对象，并设定各列宽度
    snv_table = Table(results_dict['snvList'], colWidths=[40, 60, 150, 130, 50])
    #添加表格样式
    snv_table.setStyle(TableStyle([
    ('FONTNAME',(0,0),(-1,-1),'wqy'),#字体
    ('FONTSIZE',(0,0),(-1,-1), 8),#字体大小
    ('BACKGROUND',(0,0),(-1,0), colors.lightskyblue),#设置第一行背景颜色
    ('ALIGN',(-1,0),(-2,0),'RIGHT'),#对齐
    ('VALIGN',(-1,0),(-2,0),'MIDDLE'),  #对齐
    ('LINEBEFORE',(0,0),(0,-1),0.1, colors.black),#设置表格左边线颜色为灰色，线宽为0.1
    ('GRID',(0,0),(-1,-1),0.5, colors.black),#设置表格框线为红色，线宽为0.5
    ]+snv_label_colour))
    story.append(snv_table)
    
    text = '<para autoLeading="off"><br/><font face="wqy" fontSize=8>2.基因融合</font><br/><br/></para>'
    story.append(Paragraph(text,normalStyle))

    #表格数据：用法详见reportlab-userguide.pdf中chapter 7 Table
    fusion_head = ['融合基因', '融合类型', '外显子拼接', 'Cosmic ID', '是否检出']
    results_dict['fusionList'].insert(0, fusion_head)
    #创建表格对象，并设定各列宽度
    fusion_table = Table(results_dict['fusionList'], colWidths=[40, 150, 130, 60, 50])
    #添加表格样式
    fusion_table.setStyle(TableStyle([
    ('FONTNAME',(0,0),(-1,-1),'wqy'),#字体
    ('FONTSIZE',(0,0),(-1,-1), 8),#字体大小
    ('BACKGROUND',(0,0),(-1,0), colors.lightskyblue),#设置第一行背景颜色
    ('ALIGN',(-1,0),(-2,0),'RIGHT'),#对齐
    ('VALIGN',(-1,0),(-2,0),'MIDDLE'),  #对齐
    ('LINEBEFORE',(0,0),(0,-1),0.1, colors.black),#设置表格左边线颜色为灰色，线宽为0.1
    ('TEXTCOLOR',(0,1),(-2,-1), colors.black),#设置表格内文字颜色
    ('GRID',(0,0),(-1,-1),0.5, colors.black),#设置表格框线为红色，线宽为0.5
    ]+fusion_label_colour))
    story.append(fusion_table)

    doc = SimpleDocTemplate(pdf_report)
    doc.build(story)

if __name__ == '__main__':
    json_report = sys.argv[1]
    pdf_report = sys.argv[2]
    sample_code = sys.argv[3]
    library_code = sys.argv[4]
    if len(sys.argv) == 5:
        rpt(json_report, pdf_report, sample_code, library_code)
    elif len(sys.argv) == 6:
        date = "{}年{}月{}日".format(sys.argv[5][:4],sys.argv[5][4:6],sys.argv[5][6:8])
        rpt(json_report, pdf_report, sample_code, library_code,date) 
