import subprocess
import sys

import pandas as pd
import os
from os.path import exists,abspath,dirname
import numpy as np  



def del_Unnamed(df):
    """
    Deletes all the unnamed columns

    :param df: pandas dataframe
    """
    cols_del=[c for c in df.columns if 'Unnamed' in c]
    return df.drop(cols_del,axis=1)

def genomeocoords2sections(genomecoord):   
    
    try:
        chrom=genomecoord.split(':')[0]
    except:
        raise ValueError(genomecoord)
    
    start=genomecoord.split(':')[1].split('-')[0]

    end=genomecoord.split(':')[1].split('-')[1].replace('+','').replace('-','')

    tail=genomecoord.split(':')[1].replace(start,'')
    if tail.endswith('+'):
        strand='+'
    elif tail.endswith('-'):
        strand='-'
    else:
        strand=''
#     print(tail,strand)
    return chrom,start,end,strand


def genomeocoords2bed(df, col_genomeocoord):

    bed_colns=['chromosome', 'start', 'end', 'id', 'NM', 'strand']
    df=df.dropna(subset=[col_genomeocoord])
    dbed=df.apply(lambda x: genomeocoords2sections(x[col_genomeocoord]),axis=1).apply(pd.Series)
    if len(dbed)!=0:
        dbed.columns=['chromosome', 'start', 'end','strand']
        dbed['id']=df[col_genomeocoord]
        dbed['NM']=np.nan
        return dbed.replace('',np.nan).loc[:,bed_colns].dropna(subset=['chromosome','start','end','id','strand'],how='any')
    else:
        return pd.DataFrame(columns=['chromosome', 'start', 'end','strand','id','NM'])
    
    
def fa2df(alignedfastap,ids2cols=False):
    dtmp=pd.read_csv(alignedfastap,names=["c"],keep_default_na=False)
    dtmp=dtmp.iloc[::2].reset_index(drop=True).join(dtmp.iloc[1::2].reset_index(drop=True),rsuffix='r')
    dtmp.columns=['id','sequence']
    dtmp=dtmp.set_index('id')
    dtmp.index=[i[1:] for i in dtmp.index]
    dtmp.index.name='id'
    if ids2cols:
        for i in dtmp.index:
            seqid,contig,strand,start,end=i.split('|')
            dtmp.loc[i,'seqid']=seqid
            dtmp.loc[i,'contig']=contig
            dtmp.loc[i,'strand']=strand
            dtmp.loc[i,'start']=start
            dtmp.loc[i,'end']=end
    return dtmp

def get_seq_nucleotide(genomep, dinp, bedp, fastap, dbedntmutsp, dsequencesp,flankntc):
    """
    Fetches sequences if mutation format is nucleotide

    :param cfg: configuration dict
    :param din: input data
    :returns dsequences: dataframe with sequences
    """    
  
    din = del_Unnamed(pd.read_table(dinp,keep_default_na=False) )

    if not exists(bedp):
        dbed=genomeocoords2bed(din,col_genomeocoord='genome coordinate')
        dbed['start']=dbed['start'].astype(int)-flankntc-1                          
        dbed['end']=dbed['end'].astype(int) + flankntc   
        dbed.to_csv(bedp,sep='\t',header=False, index=False)
    if not exists(fastap):
        cmd=f"bedtools getfasta -s -name -fi {genomep} -bed {bedp} -fo {fastap}"
        print(cmd)
        os.system(cmd) 
    if not exists(dbedntmutsp):
        dbedntmuts=fa2df(fastap)
        dbedntmuts.columns=['transcript: sequence']
        dbedntmuts['transcript: sequence']=dbedntmuts.apply(lambda x: x['transcript: sequence'].upper(),axis=1)
        dbedntmuts=dbedntmuts.reset_index()
        dbedntmuts['genome coordinate']=dbedntmuts.apply(lambda x : x['id'].split('(')[0] ,axis=1)
        dbedntmuts.to_csv(dbedntmutsp,sep='\t')
    else:
        dbedntmuts=pd.read_table(dbedntmutsp,keep_default_na=False)

    dsequences=del_Unnamed(pd.merge(din,dbedntmuts,on=['genome coordinate'],suffixes=('', ': dbedntmuts')))

    col_nt_wt='nucleotide wild-type' if not 'nucleotide wild-type' in dsequences else 'nucleotide wild-type: from flanking sequence'    
    col_nt_mt='nucleotide mutation' if not 'nucleotide mutation' in dsequences else 'nucleotide mutation: from flanking sequence'    
    col_cd_wt='codon: wild-type' if not 'codon: wild-type' in dsequences else 'codon: wild-type: from flanking sequence'
    col_cd_mt='codon: mutation' if not 'codon: mutation' in dsequences else 'codon: mutation: from flanking sequence'        
    dsequences[col_nt_wt]=dsequences.apply(lambda x: x['transcript: sequence'][flankntc],axis=1)        

    dsequences[col_cd_wt]=dsequences.apply(lambda x: x['transcript: sequence'][flankntc-1:flankntc+2],axis=1)
    dsequences[col_cd_mt]=dsequences.apply(lambda x: f"{x['codon: wild-type'][0]}{x['nucleotide mutation']}{x['codon: wild-type'][2]}",axis=1)
    dsequences['transcript: id']=dsequences['genome coordinate'] 

    dsequences_bedcols=genomeocoords2bed(dsequences, col_genomeocoord='genome coordinate')

    for col in dsequences_bedcols:
        dsequences[col]=dsequences_bedcols[col]
    dsequences.to_csv(dsequencesp,sep='\t')   
    
    return dsequences



def main():

    flankntc = 22
    p_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
    basedir = os.path.join(p_dir,'data/YB01/01_sequences/')

    if not exists(basedir):   
        os.makedirs(basedir)
        
    genomep = '/home/yanghe/baseEditor/beditor/temp/YB01/YB01.fna'
    dinp = '/home/yanghe/githubcode/beditor/sgRNA_primer/input_output/YB01/00_input/dinput.tsv'

    bedp = os.path.join(p_dir,'data/YB01/01_sequences/dbedntmuts.bed')
    fastap = os.path.join(p_dir,'data/YB01/01_sequences/dbedntmuts.fa')
    dbedntmutsp = os.path.join(p_dir,'data/YB01/01_sequences/dbedntmuts.tsv')
    dsequencesp = os.path.join(p_dir,'data/YB01/01_sequences/dsequences.tsv')

    dsequences = get_seq_nucleotide(genomep, dinp, bedp, fastap, dbedntmutsp, dsequencesp,flankntc)



if __name__ == '__main__':

    main()