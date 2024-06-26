import sys,os
import pandas as pd
from Bio import SeqIO 
import re
from primerDesign.src.sgRNA_utils import sgRNA_primer_util  as su
from primerDesign.src.sgRNA_module import sequencing_primer as sequencing_primer
from primerDesign.src.sgRNA_utils import sgRNA_primer_config as config     

  
#构建酶库
def create_enzyme(path='/home/yanghe/primer/editor_gene/data/enzyme.csv'):
    enzyme_df = pd.DataFrame(data=[['BsaI','GGTCTC',1,4],
                                   ['BbsI','GAAGAC',2,4],
                                   ['SapI','GCTCTTC',1,3],
                                   ['BsmBI','CGTCTC',1,4],
                                   ['BbsI','GAAGAC',2,4],
                                  ], columns=['enzyme','recognition_seq','gap_len','cut_seq_len'])
    enzyme_df.to_csv(path,index=False)  

#从质粒中搜索酶切位点
def search_enzyme_from_plasmid(enzyme_df,gb_seq_tuple,cut_num):
    strand = ''
    cut_seq = ''
    recognition_seq = ''
    gb_strand_seq,gb_re_strand_seq = gb_seq_tuple
    enzyme_in_plasmid = []
    
    for i,v in enzyme_df.iterrows():
        recognition_seq =v['recognition_seq'] 
        cut_site_cor,count = su.find_coordinate_by_pattern(recognition_seq, gb_strand_seq)
        re_cut_site_cor,re_count = su.find_coordinate_by_pattern(recognition_seq, gb_re_strand_seq)
        gap_len = v['gap_len']
        cut_seq_len = v['cut_seq_len']
        
        if count == cut_num:
            #酶切位点+ 
            temp={}
            for i,value in cut_site_cor.items(): 
                cut_start = value[1]+gap_len
                cut_end = value[1]+gap_len+cut_seq_len
                cut_coordinate = (cut_start, cut_end)
                cut_seq = gb_strand_seq[cut_start:cut_end]
                temp.update({'strand':'+','enzyme':v['enzyme'],'recognition_seq':recognition_seq,f'cut_seq_{i}':cut_seq,f'cut_coordinate_{i}':cut_coordinate})
            enzyme_in_plasmid.append(temp)
        if re_count == cut_num:
             #酶切位点-
            temp={}
            for i,value in re_cut_site_cor.items(): 
                cut_start = value[1]+gap_len
                cut_end = value[1]+gap_len+cut_seq_len
                cut_seq = gb_re_strand_seq[cut_start:cut_end]
                #取反向互补的切割序列序列
                cut_seq = su.revComp(cut_seq)
                #变换坐标
                cut_start = len(gb_re_strand_seq[cut_end:])
                cut_end = len(gb_re_strand_seq[cut_end:])+3         
                cut_coordinate = (cut_start, cut_end)
                temp.update({'strand':'-','enzyme':v['enzyme'],'recognition_seq':recognition_seq,f'cut_seq_{i}':cut_seq,f'cut_coordinate_{i}':cut_coordinate})     
            enzyme_in_plasmid.append(temp)
    return  enzyme_in_plasmid

#判断是否有可用酶,并根据用户需求选择酶------------------
def judge_is_not_availability_enzyme(enzyme_in_plasmid,enzyme_name='none'):
    enzyme_in_plasmid_df = pd.DataFrame(enzyme_in_plasmid)
    enzyme_in_plasmid_group = enzyme_in_plasmid_df.groupby('enzyme').count()

    final_enzyme_list = list(enzyme_in_plasmid_group[enzyme_in_plasmid_group['strand']==2].index)
    #默认取第一个酶
    if enzyme_name == 'none':
        final_enzyme = enzyme_in_plasmid_df[enzyme_in_plasmid_df['enzyme']==final_enzyme_list[0]]
    else:
        #酶唯一  
        final_enzyme = enzyme_in_plasmid_df[enzyme_in_plasmid_df['enzyme']==enzyme_name]
    return final_enzyme

#合并sgRNA和enzyme表
def merge_sgRNA_enzyme(final_enzyme,sgRNA_df):
    #将切割位点信息整合到sgRNA中
    sgRNA_df['enzyme']=final_enzyme.loc[1,'enzyme']
    sgRNA_df['recognition_seq'] = final_enzyme.loc[1,'recognition_seq']
    sgRNA_df['forward_joint'] = final_enzyme[final_enzyme['strand']=='-'].reset_index().loc[0,'cut_seq_0']
    sgRNA_df['forward_joint_coordinate_from_plasmid'] = str(final_enzyme[final_enzyme['strand']=='-'].reset_index().loc[0,'cut_coordinate_0'])[1:-1]
    sgRNA_df['backward_joint'] = final_enzyme[final_enzyme['strand']=='+'].reset_index().loc[0,'cut_seq_0']
    sgRNA_df['backward_joint_coordinate_from_plasmid'] = str(final_enzyme[final_enzyme['strand']=='+'].reset_index().loc[0,'cut_coordinate_0'])[1:-1]
    sgRNA_df = sgRNA_df[sgRNA_df['guide: id']!='']  
    return sgRNA_df   

#构建重组sgRNA质粒
def create_sgRNA_plasmid(sgRNA_df,gb_strand_seq): 
  
    def work(x1,x2,x3):
        arr1 = x2.split(',')   
        arr2 = x3.split(',')
        sgRNA_plasmid_seq = gb_strand_seq[:int(arr1[1].strip())] + x1 + gb_strand_seq[int(arr2[0].strip()):]
        return sgRNA_plasmid_seq
    sgRNA_df['sgRNA_plasmid_seq'] = sgRNA_df.apply(lambda x: work(x['guide+PAM sequence'],x['forward_joint_coordinate_from_plasmid'],x['backward_joint_coordinate_from_plasmid']),axis=1)
    sgRNA_df = sgRNA_df[sgRNA_df['guide: id']!='']
    return sgRNA_df


#构建重组载体的测序引物
def create_plasmid_sequencing_primer(x,y,sgRNA_df):
    sgRNA_plasmid_df = pd.DataFrame()
    for i,v in sgRNA_df.iterrows():
        sgRNA_plasmid_seq=v['sgRNA_plasmid_seq']
        sgRNA_plasmid_seq_len = len(sgRNA_plasmid_seq)

        sgRNA_len = len(v['guide+PAM sequence'])
        forward_joint_coordinate_from_plasmid_right = int(v['forward_joint_coordinate_from_plasmid'].split(',')[1])
        sgRNA_in_plasmid_coordinate = forward_joint_coordinate_from_plasmid_right,forward_joint_coordinate_from_plasmid_right+sgRNA_len
        #定义测序区域
        target_gene_seq = sgRNA_plasmid_seq[x:y]
        target_gene_up_seq = sgRNA_plasmid_seq[0:x]  
        target_gene_down_seq = sgRNA_plasmid_seq[y:sgRNA_plasmid_seq_len]

        #定义用于测序的数据结构
        dict_plasmid_seq={}
        dict_plasmid_seq['target_gene_down_seq']=target_gene_down_seq
        dict_plasmid_seq['target_gene_up_seq']=target_gene_up_seq
        dict_plasmid_seq['mute_after_target_gene_seq']=target_gene_seq     #用户指定测序序列的区域序列

        #测序质粒的id
        result = {'guide: id':v['guide: id']}
        result.update(sequencing_primer.design_sequencing_primers(v['guide: id'], dict_plasmid_seq)[0])
        sgRNA_plasmid_df = sgRNA_plasmid_df.append(pd.DataFrame([result]))
    sgRNA_plasmid_df = su.make_id_in_first_columns(sgRNA_plasmid_df,id='guide: id',columns=sgRNA_plasmid_df.columns.tolist())
    
    return sgRNA_plasmid_df

#取出基因组突变位点中的前后bp
def extract_forward_backward_from_genome(sgRNA_df, genome_fna_path, mutation_info):
    # mutation_info.drop(columns=['nucleotide mutation'], inplace=True) 
    sgRNA_df["genome coordinate"] = sgRNA_df["guide: id"].apply(lambda x: x.split("|")[0] )

    sgRNA_df = pd.merge(sgRNA_df,mutation_info,on='genome coordinate',how='inner')  
    sgRNA_df["forward_seq_bp"] = 300  
    sgRNA_df["backward_seq_bp"] = 300    

    def work(x1,x2,x3,x4):
        
        arr = x1.split(':')
        id = arr[0]
        position = int(arr[1].split('-')[0])   
        strand = arr[1][-1]
        for record in SeqIO.parse(genome_fna_path, "fasta"):
            if record.id == id:
                if strand == '+':   
                    mute = x2  
                elif strand =='-':
                    mute = su.revComp(x2) 
            
                length = len(str(record.seq))
                x3 = int(x3)    
                x4 = int(x4)
                if x3 < 200:
                    temp_seq =  str(record.seq[-(200-x3):])
                    target_seq_up = temp_seq + str(record.seq)[:x3]
                elif length - x4 < 200:
                    temp_seq =  str(record.seq[:200 - (length - x4)])
                    target_seq_down = str(record.seq)[x4:]+temp_seq
                else:
                    target_seq_up = str(record.seq)[x3-200:x3]
                    target_seq_down = str(record.seq)[x4:x4+200]
                target_seq = str(record.seq)[x3:x4+1]
            return f'{target_seq_up},{target_seq},{target_seq_down}'   
    sgRNA_df['forward_seq_backward']=sgRNA_df.apply(lambda x:work(x['genome coordinate'],x['nucleotide mutation'],x['forward_seq_bp'],x['backward_seq_bp']),axis=1)  

   
    return sgRNA_df

  
#创建编辑序列的测序引物
def create_editor_sequencing_primer(sgRNA_df):   
    editor_sequencing = pd.DataFrame()
    for i,v in sgRNA_df.iterrows():
        id = v['genome coordinate']
        arr = v['forward_seq_backward'].split(',')
        left_seq = arr[0]
        safe_area = arr[1]
        right_seq =  arr[2]  

        #定义用于测序的数据结构
        dict_plasmid_seq={}
        dict_plasmid_seq['target_gene_down_seq']=right_seq
        dict_plasmid_seq['target_gene_up_seq']=left_seq
        dict_plasmid_seq['mute_after_target_gene_seq']=safe_area     #用户指定测序序列的区域序列

         #测序质粒的id
        result = {'guide: id':v['guide: id']}
        success,failtrue = sequencing_primer.design_sequencing_primers(v['guide: id'],dict_plasmid_seq)

        result.update(success)
      
        editor_sequencing = editor_sequencing.append(pd.DataFrame([result]))
       
    editor_sequencing = su.make_id_in_first_columns(df=editor_sequencing,id='guide: id',columns=editor_sequencing.columns.tolist())
    return editor_sequencing  


#判断引物是与否
# def judge_primer_is_or_not(dict_res,primer_suc,primer_fail,primer_name,type='LEFT'):
#     if len(list(dict_res.keys())) < 10:
#         primer_fail[primer_name] = dict_res   
#     else:
#         primer_suc[primer_name] = dict_res[f'PRIMER_{type}_0_SEQUENCE'] 

##构建sgRNA的oligo
def create_sgRNA_oligo(sgRNA_df):
    sgRNA_df['sgRNA_oligo+']=sgRNA_df['forward_joint']+sgRNA_df['guide+PAM sequence']
    sgRNA_df['sgRNA_oligo-']=sgRNA_df.apply(lambda x: x['backward_joint'] + su.revComp(x['guide+PAM sequence']),axis=1)
    sgRNA_oligo_df = sgRNA_df[['guide: id','sgRNA_oligo+','sgRNA_oligo-']]  
    return sgRNA_oligo_df 

#添加
   

def main():
    #加载gb文件读取序列信息
    gb = SeqIO.read('/home/yanghe/githubcode/beditor/sgRNA_primer/temp/alsD-28a1.gb', "genbank")   
    gb_strand_seq = str(gb.seq)
    gb_re_strand_seq = su.revComp(gb_strand_seq)
    gb_seq_tuple = (gb_strand_seq,gb_re_strand_seq)
    
    #加载酶信息
    enzyme_df = pd.read_csv('/home/yanghe/githubcode/beditor/sgRNA_primer/input_output/enzyme.csv') 
    #搜索可用酶
    enzyme_in_plasmid = search_enzyme_from_plasmid(enzyme_df,gb_seq_tuple,cut_num=1)
    #判断是否有可用酶,并根据用户需求选择酶
    final_enzyme = judge_is_not_availability_enzyme(enzyme_in_plasmid,enzyme_name='none')
    assert len(final_enzyme)!=0 
      #将酶与sgRNA_df两个表整合
    sgRNA_df = merge_sgRNA_enzyme(final_enzyme,sgRNA_df)
    #构建sgRNA的oligo
    sgRNA_oligo_df = create_sgRNA_oligo(sgRNA_df)
    #构建sgRNA重组质粒
    sgRNA_plasmid_df=create_sgRNA_plasmid(sgRNA_df,gb_strand_seq)
    #构建重组载体的测序引物
    #用户自定义扩展sgRNA的坐标，0<x<2261,len(plasmid_seg)>y>2284
    x = 2240 
    y = 5000
    plasmid_sequencing_primer_df = create_plasmid_sequencing_primer(x,y,sgRNA_plasmid_df,sequencing_primer)
    
    #构建碱基编辑的测序引物
        #取出基因组上突变位点的前后600bp
    genome_fna_path="/home/yanghe/githubcode/beditor/super_beditor/pub/GCF_000204275.1/fasta/saccharomyces_cerevisiae/dna/genome.fna" 
    sgRNA_plasmid_600bp_df = extract_forward_backward_bp_from_genome(sgRNA_plasmid_df,genome_fna_path) 

    editor_sequencing_df = create_editor_sequencing_primer(sgRNA_plasmid_600bp_df)
    
    #输出
    sheetname2df = {'sgRNA_info':sgRNA_df,'sgRNA_oligo':sgRNA_oligo_df,'editor_sequencing':editor_sequencing_df,'plasmid_sequencing':plasmid_sequencing_primer_df}
    datap = "/home/yanghe/githubcode/beditor/sgRNA_primer/input_output/superBeditor_sgRNA_primer_output.xlsx"
    su.to_excel(sheetname2df,datap,)