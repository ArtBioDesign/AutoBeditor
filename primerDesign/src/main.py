from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
import os,sys 
import pandas as pd
from Bio import SeqIO 
from primerDesign.src.sgRNA_utils import sgRNA_primer_util as su
from primerDesign.src.sgRNA_module import desgin_sgRNA_primer as sr


def oligo_design(gb,enzyme_df,sgRNA_df, genome_path, mutation_info, workdir):
    
    #1 加载gb文件读取序列信息
    gb_strand_seq = str(gb.seq)
    gb_re_strand_seq = su.revComp(gb_strand_seq)
    gb_seq_tuple = (gb_strand_seq, gb_re_strand_seq)

    
    #2 搜索可用酶
    enzyme_in_plasmid = sr.search_enzyme_from_plasmid(enzyme_df, gb_seq_tuple, cut_num=1)
    #判断是否有可用酶,并根据用户需求选择酶
    final_enzyme = sr.judge_is_not_availability_enzyme(enzyme_in_plasmid, enzyme_name='none')
    assert len(final_enzyme)!=0 

    #3. 将酶与sgRNA_df两个表整合
    sgRNA_df = sr.merge_sgRNA_enzyme(final_enzyme,sgRNA_df)

    #4 构建sgRNA的oligo
    sgRNA_oligo_df = sr.create_sgRNA_oligo(sgRNA_df)
    #5 构建sgRNA重组质粒
    sgRNA_plasmid_df=sr.create_sgRNA_plasmid(sgRNA_df,gb_strand_seq)

    #6 构建重组载体的测序引物
    #用户自定义扩展sgRNA的坐标，0<x<2261,len(plasmid_seg)>y>2284
    x = 2240 
    y = 5000
    plasmid_sequencing_primer_df = sr.create_plasmid_sequencing_primer(x,y,sgRNA_plasmid_df)
    
    #7 构建碱基编辑的测序引物
        #取出基因组上突变位点的前后600bp
        # extract_forward_backward_from_genome(sgRNA_df, genome_fna_path, mutation_info)
    sgRNA_plasmid_600bp_df = sr.extract_forward_backward_from_genome(sgRNA_plasmid_df, genome_path, mutation_info) 
    editor_sequencing_df = sr.create_editor_sequencing_primer(sgRNA_plasmid_600bp_df)
    
    #8 输出
    sheetname2df = {'sgRNA_info':sgRNA_df,'sgRNA_oligo':sgRNA_oligo_df,'editor_sequencing':editor_sequencing_df,'plasmid_sequencing':plasmid_sequencing_primer_df}
    datap = os.path.join(workdir, "sgRNA_primer_output.xlsx")
    su.to_excel(sheetname2df,datap)


def run(sgRNA_path, gb_path, enzy_path, genome_path, mutation_info_path, workdir):

    #读取sgRNA
    sgRNA_df = su.del_Unnamed(pd.read_csv(sgRNA_path ,sep='\t',keep_default_na=False))
    mutation_info = su.del_Unnamed(pd.read_csv(mutation_info_path, sep='\t',keep_default_na=False))

    
    #设计相关的oligo、primer信息
    gb = SeqIO.read(gb_path, "genbank")
    enzyme_df = su.del_Unnamed( pd.read_csv(enzy_path) )
    oligo_design(gb,enzyme_df,sgRNA_df, genome_path, mutation_info, workdir)


def main():
    cur = "/hpcfs/fhome/yangchh/editSeqDesign/self/AutoBeditor"
    sgRNA_path = os.path.join(cur, "data/output/yb01/04_offtargets/dofftargets.tsv")
    gb_path  = os.path.join(cur, "data/input/alsD-28a1.gb")
    enzy_path = os.path.join(cur, 'data/input/enzyme.csv')
    genome_path = os.path.join(cur, "data/input/YB01.fna")
    workdir = os.path.join(cur, "data/output/yb01")
    mutation_info_path = os.path.join(cur, "data/input/YB01_beditor_input1.tsv")
    run(sgRNA_path, gb_path, enzy_path, genome_path, mutation_info_path, workdir)
    

if __name__ == '__main__':
     
     main()
