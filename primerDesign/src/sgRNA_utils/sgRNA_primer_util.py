import re
import pandas as pd
import primer3
from primerDesign.src.sgRNA_utils.sgRNA_primer_config import config   
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
import os

def del_Unnamed(df):
    """
    Deletes all the unnamed columns

    :param df: pandas dataframe
    """
    cols_del=[c for c in df.columns if 'Unnamed' in c]
    return df.drop(cols_del,axis=1)

def dfswapcols(df,cols):                #交换列
    df[f"_{cols[0]}"]=df[cols[0]].copy()
    df[cols[0]]=df[cols[1]].copy()
    df[cols[1]]=df[f"_{cols[0]}"].copy()
    df=df.drop([f"_{cols[0]}"],axis=1)
    return df

#取互补链
def revComp(seq):
    complementSeq=seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
    revcompSeq = complementSeq[::-1]
    return revcompSeq

#发现字符串的位置
def find_coordinate_by_pattern(pattern,seq):
    sub_seq_cor = {}
    i = 0
    for m in re.finditer(pattern, seq):
        sub_seq_cor.update({f'{i}':(m.start(), m.end())})
        i = i + 1
    return sub_seq_cor,i

#使id位于第一列
def make_id_in_first_columns(df,id,columns):
    assert id in columns
    first_columns_list=df[id].tolist()
    df.drop(columns = [id],inplace=True)
    df.insert(loc=0,column =id ,value = first_columns_list)
    return df

def read_excel(p,sheet_name=None,):
    xl = pd.ExcelFile(p)
    xl.sheet_names  # see all sheet names
    if sheet_name is None:
        sheet_name=input(', '.join(xl.sheet_names))
    return xl.parse(sheet_name) 

def to_excel(sheetname2df,datap,):
    
    if not exists( dirname(datap) ):
        os.makedirs( dirname(datap) )

    writer = pd.ExcelWriter(datap)
    for sn in sheetname2df:
        sheetname2df[sn].to_excel(writer,sn)
    writer.save() 


#换名字
def replace_primer3Name_to_peopleReadName(df,type='up'):
    names = df.columns.to_list()
    df =df.rename(columns={
                            names[2]:f"primer_{type}_f_seq_(5'-3')",
                            names[3]:f"primer_{type}_r_seq_(5'-3')",
                            names[4]:f"primer_{type}_f_Tm",
                            names[5]:f"primer_{type}_r_Tm",
                            names[1]:f"{type}_product_sequence_length",
                            names[0]:f"{type}_product_sequence"
                        } )
    return df  

#设计引物
def primer_design(seqId,
                  seqTemplate,
                  stype,
                  mute_type='single',
                  global_args=config.GLOBAL_ARGS):

    if mute_type == 'single':
        config.GLOBAL_ARGS.update(config.S_GLOBAL_ARGS)
        global_args = config.GLOBAL_ARGS 
    elif mute_type =='double':
        config.GLOBAL_ARGS.update(config.D_GLOBAL_ARGS)
        global_args = config.GLOBAL_ARGS
    elif mute_type == 'sequencing': 
        config.Q_GLOBAL_ARGS.update(config.Q_ARGS)
        global_args =config.Q_GLOBAL_ARGS 

        
    #序列参数  
    seqlength = len(seqTemplate)   
    seq_args = {
                'SEQUENCE_ID': seqId,
                'SEQUENCE_TEMPLATE': seqTemplate,
                'SEQUENCE_FORCE_LEFT_START':-1,
                'SEQUENCE_FORCE_RIGHT_START': -1
                # 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [[0,0,0,0]]
        }
    #选择类型，设定序列参数
    if mute_type == 'single': #单点突变
        if stype == "left":   #target上游引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[0,100,-1,-1]]
            size_range = [seqlength-66,seqlength-7]
        elif stype == "right":   #target下游引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[-1,-1,seqlength-40,40]]
            seq_args['SEQUENCE_FORCE_LEFT_START'] = 0
            size_range = [seqlength-2,seqlength]
        primer_num = 2    
    elif mute_type == 'double':  #双点突变                                   
        if stype == "target":   #target引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[0,35,-1,-1]]
        elif stype == "plasmid":   #plasmid引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[-1,-1,seqlength-36,36]]     
        size_range = [seqlength,seqlength]
        primer_num = 2
    elif mute_type == 'sequencing':  #测序引物
        seq_args['SEQUENCE_ID'] = seqId
        seq_args['SEQUENCE_TEMPLATE'] = seqTemplate  
        size_range = [25,seqlength]
        primer_num = 1

    #设定全局参数   
    global_args['PRIMER_PRODUCT_SIZE_RANGE']=size_range
    global_args['PRIMER_NUM_RETURN']= primer_num
  
    #调用工具
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args)    
    return primer3_result

#输出引物设计成功的
def result_output_success_df(plasmid_name,primers_dict,type='down'):
    all_df = pd.DataFrame()
    
    for key in primers_dict:
        df =pd.DataFrame(primers_dict[key])    
        df['id'] = key   
        all_df=all_df.append(df)  
    #换名子
    all_df = replace_primer3Name_to_peopleReadName(all_df,type)  
    columns_name = all_df.columns.to_list()
    # #列的顺序
    all_df = all_df[[columns_name[-1],columns_name[2], columns_name[3],columns_name[4],columns_name[5],columns_name[0],columns_name[1]]]
    all_df['template']=plasmid_name
    return all_df  

#输出设计引物失败  
def result_output_failtrue_df(plasmid_name,primers_dict_failture):

    df = pd.DataFrame()
    for k,v in primers_dict_failture.items():
        temp = pd.DataFrame([v])
        temp.insert(0,'id',k)
        df = df.append(temp)
    df['template']=plasmid_name
    return df  