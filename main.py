
from beditor import main as sgRNA_main
from primerDesign.src import main as primer_main
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
import os,sys 
import pandas as pd

def main():

    #1.上传数据文件
    cur = "/hpcfs/fhome/yangchh/editSeqDesign/self/AutoBeditor"
    genome_fna_path = os.path.join(cur, "data/input/YB01.fna")
    genome_gff_path = os.path.join(cur, "data/input/YB01_v2.gff")
    gb = os.path.join(cur, "data/input/alsD-28a1.gb")
    dinp = os.path.join(cur, "data/input/YB01_beditor_input1.tsv") 
    enzy =  os.path.join(cur, 'data/input/enzyme.csv')
   
    #2.碱基编辑设计sgRNA
    workdir = os.path.join(cur,'data/output')  
    if not exists(workdir):
        os.makedirs( workdir ) 

    host = 'Yb01'
    param={ 
        "genome_release" : 'GCF001',
        "host" : host,
        "be_names" : ['Target-AID'],
        "pams" : ['NG']
    } 
    sgRNA_main.run(dinp, genome_fna_path, genome_gff_path, workdir, param)
  
    #5.读取sgRNA
    sgRNA_path = os.path.join(cur, f"data/output/{host}/04_offtargets/dofftargets.tsv")

    # 6.设计相关的oligo、primer信息
    workdir = os.path.join(cur,f'data/output/{host}')
    primer_main.run(sgRNA_path, gb, enzy, genome_fna_path, dinp, workdir)
    
      
if __name__ == '__main__':
     
     main()
