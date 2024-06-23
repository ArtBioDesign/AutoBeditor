import os
import pandas as pd
import logging
from glob import glob
from os import makedirs,stat
import pysam
import numpy as np
import re
from os.path import join, basename, dirname, abspath, exists

def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

def calc_cfd(wt,sg,pam):
    mm_scores = {'rA:dA,1': 1.0,
 'rA:dA,10': 0.882352941,
 'rA:dA,11': 0.307692308,
 'rA:dA,12': 0.333333333,
 'rA:dA,13': 0.3,
 'rA:dA,14': 0.533333333,
 'rA:dA,15': 0.2,
 'rA:dA,16': 0.0,
 'rA:dA,17': 0.133333333,
 'rA:dA,18': 0.5,
 'rA:dA,19': 0.538461538,
 'rA:dA,2': 0.727272727,
 'rA:dA,20': 0.6,
 'rA:dA,3': 0.705882353,
 'rA:dA,4': 0.636363636,
 'rA:dA,5': 0.363636364,
 'rA:dA,6': 0.7142857140000001,
 'rA:dA,7': 0.4375,
 'rA:dA,8': 0.428571429,
 'rA:dA,9': 0.6,
 'rA:dC,1': 1.0,
 'rA:dC,10': 0.5555555560000001,
 'rA:dC,11': 0.65,
 'rA:dC,12': 0.7222222220000001,
 'rA:dC,13': 0.6521739129999999,
 'rA:dC,14': 0.46666666700000003,
 'rA:dC,15': 0.65,
 'rA:dC,16': 0.192307692,
 'rA:dC,17': 0.176470588,
 'rA:dC,18': 0.4,
 'rA:dC,19': 0.375,
 'rA:dC,2': 0.8,
 'rA:dC,20': 0.764705882,
 'rA:dC,3': 0.611111111,
 'rA:dC,4': 0.625,
 'rA:dC,5': 0.72,
 'rA:dC,6': 0.7142857140000001,
 'rA:dC,7': 0.705882353,
 'rA:dC,8': 0.7333333329999999,
 'rA:dC,9': 0.666666667,
 'rA:dG,1': 0.857142857,
 'rA:dG,10': 0.333333333,
 'rA:dG,11': 0.4,
 'rA:dG,12': 0.263157895,
 'rA:dG,13': 0.21052631600000002,
 'rA:dG,14': 0.214285714,
 'rA:dG,15': 0.272727273,
 'rA:dG,16': 0.0,
 'rA:dG,17': 0.176470588,
 'rA:dG,18': 0.19047619,
 'rA:dG,19': 0.20689655199999998,
 'rA:dG,2': 0.7857142859999999,
 'rA:dG,20': 0.22727272699999998,
 'rA:dG,3': 0.428571429,
 'rA:dG,4': 0.352941176,
 'rA:dG,5': 0.5,
 'rA:dG,6': 0.454545455,
 'rA:dG,7': 0.4375,
 'rA:dG,8': 0.428571429,
 'rA:dG,9': 0.571428571,
 'rC:dA,1': 1.0,
 'rC:dA,10': 0.9411764709999999,
 'rC:dA,11': 0.307692308,
 'rC:dA,12': 0.538461538,
 'rC:dA,13': 0.7,
 'rC:dA,14': 0.7333333329999999,
 'rC:dA,15': 0.066666667,
 'rC:dA,16': 0.307692308,
 'rC:dA,17': 0.46666666700000003,
 'rC:dA,18': 0.642857143,
 'rC:dA,19': 0.46153846200000004,
 'rC:dA,2': 0.9090909090000001,
 'rC:dA,20': 0.3,
 'rC:dA,3': 0.6875,
 'rC:dA,4': 0.8,
 'rC:dA,5': 0.636363636,
 'rC:dA,6': 0.9285714290000001,
 'rC:dA,7': 0.8125,
 'rC:dA,8': 0.875,
 'rC:dA,9': 0.875,
 'rC:dC,1': 0.913043478,
 'rC:dC,10': 0.38888888899999996,
 'rC:dC,11': 0.25,
 'rC:dC,12': 0.444444444,
 'rC:dC,13': 0.13636363599999998,
 'rC:dC,14': 0.0,
 'rC:dC,15': 0.05,
 'rC:dC,16': 0.153846154,
 'rC:dC,17': 0.058823529000000006,
 'rC:dC,18': 0.133333333,
 'rC:dC,19': 0.125,
 'rC:dC,2': 0.695652174,
 'rC:dC,20': 0.058823529000000006,
 'rC:dC,3': 0.5,
 'rC:dC,4': 0.5,
 'rC:dC,5': 0.6,
 'rC:dC,6': 0.5,
 'rC:dC,7': 0.470588235,
 'rC:dC,8': 0.642857143,
 'rC:dC,9': 0.6190476189999999,
 'rC:dT,1': 1.0,
 'rC:dT,10': 0.8666666670000001,
 'rC:dT,11': 0.75,
 'rC:dT,12': 0.7142857140000001,
 'rC:dT,13': 0.384615385,
 'rC:dT,14': 0.35,
 'rC:dT,15': 0.222222222,
 'rC:dT,16': 1.0,
 'rC:dT,17': 0.46666666700000003,
 'rC:dT,18': 0.538461538,
 'rC:dT,19': 0.428571429,
 'rC:dT,2': 0.727272727,
 'rC:dT,20': 0.5,
 'rC:dT,3': 0.8666666670000001,
 'rC:dT,4': 0.842105263,
 'rC:dT,5': 0.571428571,
 'rC:dT,6': 0.9285714290000001,
 'rC:dT,7': 0.75,
 'rC:dT,8': 0.65,
 'rC:dT,9': 0.857142857,
 'rG:dA,1': 1.0,
 'rG:dA,10': 0.8125,
 'rG:dA,11': 0.384615385,
 'rG:dA,12': 0.384615385,
 'rG:dA,13': 0.3,
 'rG:dA,14': 0.26666666699999997,
 'rG:dA,15': 0.14285714300000002,
 'rG:dA,16': 0.0,
 'rG:dA,17': 0.25,
 'rG:dA,18': 0.666666667,
 'rG:dA,19': 0.666666667,
 'rG:dA,2': 0.636363636,
 'rG:dA,20': 0.7,
 'rG:dA,3': 0.5,
 'rG:dA,4': 0.363636364,
 'rG:dA,5': 0.3,
 'rG:dA,6': 0.666666667,
 'rG:dA,7': 0.571428571,
 'rG:dA,8': 0.625,
 'rG:dA,9': 0.533333333,
 'rG:dG,1': 0.7142857140000001,
 'rG:dG,10': 0.4,
 'rG:dG,11': 0.428571429,
 'rG:dG,12': 0.529411765,
 'rG:dG,13': 0.42105263200000004,
 'rG:dG,14': 0.428571429,
 'rG:dG,15': 0.272727273,
 'rG:dG,16': 0.0,
 'rG:dG,17': 0.235294118,
 'rG:dG,18': 0.47619047600000003,
 'rG:dG,19': 0.448275862,
 'rG:dG,2': 0.692307692,
 'rG:dG,20': 0.428571429,
 'rG:dG,3': 0.384615385,
 'rG:dG,4': 0.529411765,
 'rG:dG,5': 0.7857142859999999,
 'rG:dG,6': 0.681818182,
 'rG:dG,7': 0.6875,
 'rG:dG,8': 0.615384615,
 'rG:dG,9': 0.538461538,
 'rG:dT,1': 0.9,
 'rG:dT,10': 0.933333333,
 'rG:dT,11': 1.0,
 'rG:dT,12': 0.933333333,
 'rG:dT,13': 0.923076923,
 'rG:dT,14': 0.75,
 'rG:dT,15': 0.9411764709999999,
 'rG:dT,16': 1.0,
 'rG:dT,17': 0.933333333,
 'rG:dT,18': 0.692307692,
 'rG:dT,19': 0.7142857140000001,
 'rG:dT,2': 0.846153846,
 'rG:dT,20': 0.9375,
 'rG:dT,3': 0.75,
 'rG:dT,4': 0.9,
 'rG:dT,5': 0.8666666670000001,
 'rG:dT,6': 1.0,
 'rG:dT,7': 1.0,
 'rG:dT,8': 1.0,
 'rG:dT,9': 0.642857143,
 'rU:dC,1': 0.956521739,
 'rU:dC,10': 0.5,
 'rU:dC,11': 0.4,
 'rU:dC,12': 0.5,
 'rU:dC,13': 0.260869565,
 'rU:dC,14': 0.0,
 'rU:dC,15': 0.05,
 'rU:dC,16': 0.346153846,
 'rU:dC,17': 0.117647059,
 'rU:dC,18': 0.333333333,
 'rU:dC,19': 0.25,
 'rU:dC,2': 0.84,
 'rU:dC,20': 0.176470588,
 'rU:dC,3': 0.5,
 'rU:dC,4': 0.625,
 'rU:dC,5': 0.64,
 'rU:dC,6': 0.571428571,
 'rU:dC,7': 0.588235294,
 'rU:dC,8': 0.7333333329999999,
 'rU:dC,9': 0.6190476189999999,
 'rU:dG,1': 0.857142857,
 'rU:dG,10': 0.533333333,
 'rU:dG,11': 0.666666667,
 'rU:dG,12': 0.947368421,
 'rU:dG,13': 0.7894736840000001,
 'rU:dG,14': 0.28571428600000004,
 'rU:dG,15': 0.272727273,
 'rU:dG,16': 0.666666667,
 'rU:dG,17': 0.705882353,
 'rU:dG,18': 0.428571429,
 'rU:dG,19': 0.275862069,
 'rU:dG,2': 0.857142857,
 'rU:dG,20': 0.090909091,
 'rU:dG,3': 0.428571429,
 'rU:dG,4': 0.647058824,
 'rU:dG,5': 1.0,
 'rU:dG,6': 0.9090909090000001,
 'rU:dG,7': 0.6875,
 'rU:dG,8': 1.0,
 'rU:dG,9': 0.923076923,
 'rU:dT,1': 1.0,
 'rU:dT,10': 0.857142857,
 'rU:dT,11': 0.75,
 'rU:dT,12': 0.8,
 'rU:dT,13': 0.692307692,
 'rU:dT,14': 0.6190476189999999,
 'rU:dT,15': 0.578947368,
 'rU:dT,16': 0.9090909090000001,
 'rU:dT,17': 0.533333333,
 'rU:dT,18': 0.666666667,
 'rU:dT,19': 0.28571428600000004,
 'rU:dT,2': 0.846153846,
 'rU:dT,20': 0.5625,
 'rU:dT,3': 0.7142857140000001,
 'rU:dT,4': 0.47619047600000003,
 'rU:dT,5': 0.5,
 'rU:dT,6': 0.8666666670000001,
 'rU:dT,7': 0.875,
 'rU:dT,8': 0.8,
 'rU:dT,9': 0.9285714290000001}
    pam_scores = {'AA': 0.0,
 'AC': 0.0,
 'AG': 0.25925925899999996,
 'AT': 0.0,
 'CA': 0.0,
 'CC': 0.0,
 'CG': 0.107142857,
 'CT': 0.0,
 'GA': 0.06944444400000001,
 'GC': 0.022222222000000003,
 'GG': 1.0,
 'GT': 0.016129031999999998,
 'TA': 0.0,
 'TC': 0.0,
 'TG': 0.038961038999999996,
 'TT': 0.0}
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

stepi2cols={
1: ['aminoacid: position', 'gene: id', 'gene: name', 'protein: id', 'transcript: id', 'transcript: sequence', 'aminoacid: wild-type', 'codon: wild-type', 'contig', 'strand', 'start', 'end', 'codon start', 'codon end'],
2:  ['amino acid',
 'amino acid mutation',
 'codon',
 'codon mutation',
 'codon mutation usage Fraction',
 'codon mutation usage Frequency',
 'method',
 'mutation on strand',
 'nucleotide',
 'nucleotide mutation',
 'nucleotide mutation: count',
 'position of mutation in codon'],    
3: ['PAM','Description','guide length','original','original position','PAM position','rPAM','reverse complement','strand','is a reverse complement','guide+PAM sequence','guide sequence','PAM sequence','position of PAM ini','position of PAM end','distance of codon from PAM','codon','guide sequence length','index','transcript: id','aminoacid: position','aminoacid: wild-type','codon: wild-type','id','amino acid mutation','method','amino acid','position of mutation in codon','codon mutation','nucleotide','nucleotide mutation','mutation on strand','codon mutation usage Fraction','codon mutation usage Frequency','nucleotide mutation: count','distance of mutation from PAM: minimum','distance of mutation from PAM: maximum','distance of codon start from PAM: minimum','distance of codon start from PAM: maximum','activity sequence','distance of mutation from PAM','strategy','guide: id','guide+PAM length'],
4: ['guide: id','id','guide+PAM sequence','gene names','gene ids','transcript ids','types','protein ids','exon ids','beditor score','CFD score','beditor score (log10)','alternate alignments count'],    
5:['transcript: id',
 'aminoacid: position',
 'amino acid mutation',
 'aminoacid: wild-type',
 'guide: id',
 'guide+PAM sequence',
 'beditor score',
 'alternate alignments count',
 'CFD score'],
}

def get_beditorscore_per_guide(guide_seq, strategy,
                               align_seqs_scores,
                              dBEs,
                              penalty_activity_window=0.5,
                               test=False,
                              ):
    """
    Calculates beditor score per guide.
    
    :param guide_seq: guide seqeunce 23nts
    :param strategy: strategy string eg. ABE;+;@-14;ACT:GCT;T:A;
    :param align_seqs_scores: list of beditor scores per alignments for all the alignments between guide and genomic DNA
    :param penalty_activity_window: if editable base is not in activity window, penalty_activity_window=0.5
    :returns: beditor score per guide.
    """
    
    #create pos_muts for back-compatibility
    pos_muts=dBEs.loc[:,['method','distance of mutation from PAM: minimum',
     'distance of mutation from PAM: maximum',
     'distance of codon start from PAM: minimum',
     'distance of codon start from PAM: maximum']].drop_duplicates().set_index('method')

    pos_mut=int(strategy.split(';')[2].replace('@',''))
    method=strategy.split(';')[0]
    penalty_activity_window=1 if (pos_muts.loc[method,'distance of mutation from PAM: minimum']<=pos_mut<=pos_muts.loc[method,'distance of mutation from PAM: maximum']) else penalty_activity_window
    penalty_align_seqs_scores=np.prod(align_seqs_scores)
    if test:
        print(list(align_seqs_scores))
        print(penalty_align_seqs_scores)
    return penalty_activity_window*penalty_align_seqs_scores

def lambda2cols(df,lambdaf,in_coln,to_colns):         #apply函数的助手
    df_=df.apply(lambda x: lambdaf(x[in_coln]),
                 axis=1).apply(pd.Series)
    df_.columns=to_colns
    df=df.join(df_)        
    return df

def gffatributes2ids(s):
    print(s)
    """
    Deconvolutes ids from `attributes` column in GFF3 to seprate columns.
    :param s: attribute string.
    :returns: tuple of ids
    """
    Name,gene_id,transcript_id,protein_id,exon_id=np.nan,np.nan,np.nan,np.nan,np.nan
    if not pd.isnull(s):
        if '=' in s:
            d=dict([i.split('=') for i in s.split(';')])
            if 'Parent' in d:
                d[d['Parent'].split(':')[0]+'_id']= d['Parent']                                            #d['Parent'].split(':')[1]
            Name,gene_id,transcript_id,protein_id,exon_id=np.nan,np.nan,np.nan,np.nan,np.nan
            if 'Name' in d:    
                Name=d['Name']
            if 'gene_id' in d:    
                gene_id=d['gene_id']                                            #基因的id
            if 'transcript_id' in d:                                    
                transcript_id=d['transcript_id']                             #transcript_id ---- transl_table
            if 'protein_id' in d:    
                protein_id=d['protein_id']
            if 'exon_id' in d:    
                exon_id=d['exon_id']                                     #ensembl中的exon_id ncbi中exon_number
    return Name,gene_id,transcript_id,protein_id,exon_id

def align(s1,s2,test=False,
         psm=2,pmm=0.5,pgo=-3,pge=-1):
    """
    Creates pairwise local alignment between seqeunces.
    Get the visualization and alignment scores.
    :param s1: seqeunce 1
    :param s2: seqeunce 2    
    
    REF: http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
    The match parameters are:

    CODE  DESCRIPTION
    x     No parameters. Identical characters have score of 1, otherwise 0.
    m     A match score is the score of identical chars, otherwise mismatch
          score.
    d     A dictionary returns the score of any pair of characters.
    c     A callback function returns scores.
    The gap penalty parameters are:

    CODE  DESCRIPTION
    x     No gap penalties.
    s     Same open and extend gap penalties for both sequences.
    d     The sequences have different open and extend gap penalties.
    c     A callback function returns the gap penalties.    
    """
    import operator
    from Bio import pairwise2
    if any([p is None for p in [psm,pmm,pgo,pge]]):
        alignments = pairwise2.align.localxx(s1.upper(),s2.upper())
    else:
        alignments = pairwise2.align.localms(s1.upper(),s2.upper(),psm,pmm,pgo,pge)
    if test:
        print(alignments)
    alignsymb=np.nan
    score=np.nan
    sorted_alignments = sorted(alignments, key=operator.itemgetter(2))
    for a in alignments:
        alignstr=pairwise2.format_alignment(*a)
        alignsymb=alignstr.split('\n')[1]
        score=a[2]
        if test:
            print(alignstr)
        break
    return alignsymb.replace(' ','-'),score

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

def str2num(x):
    """
    This extracts numbers from strings. eg. 114 from M114R.

    :param x: string
    """
    return int(''.join(ele for ele in x if ele.isdigit()))

def str2numorstr(x,method=int):
    """
    This extracts numbers from strings. eg. 114 from M114R.

    :param x: string
    """
    try:
        x=method(x)
        return x
    except:
        return x

def del_Unnamed(df):
    """
    Deletes all the unnamed columns

    :param df: pandas dataframe
    """
    cols_del=[c for c in df.columns if 'Unnamed' in c]
    return df.drop(cols_del,axis=1)
def saveemptytable(step,doutp=None):
   
    dout=pd.DataFrame(columns=stepi2cols[step])
    logging.warning(f"saved enmpty table {doutp}")
    if not doutp is None:
        dout.to_csv(doutp,sep='\t')
    else:
        return dout 

def df2info(df,col_searches=None):
    if len(df.columns)>5:
        print('**COLS**: ',df.columns.tolist())
    print('**HEAD**: ',df.loc[:,df.columns[:5]].head())
    print('**SHAPE**: ',df.shape)
    if not col_searches is None:
        cols_searched=[c2 for c1 in col_searches for c2 in df if c1 in c2]
        print('**SEARCHEDCOLS**:\n',cols_searched)
        print('**HEAD**: ',df.loc[:,cols_searched].head())    
    
def set_index(data,col_index):
    """
    Sets the index if the index is not present

    :param data: pandas table 
    :param col_index: column name which will be assigned as a index
    """
    if col_index in data:
        data=data.reset_index().set_index(col_index)
        if 'index' in data:
            del data['index']
        return data
    elif data.index.name==col_index:
        return data
    else:
        logging.error("something's wrong with the df")
        df2info(data)
                
def dguides2guidessam(datatmpd,genomep,dguides):
    """
    Aligns guides to genome and gets SAM file
    step#1

    :param cfg: configuration dict
    :param dguides: dataframe of guides
    """
    
    dguides=set_index(dguides,'guide: id')
    guidels=dguides.loc[:,'guide+PAM length'].unique()

    for guidel in guidels:
        logging.debug(f"now aligning guides of length {guidel}")
        guidesfap = os.path.join( datatmpd, f'01_guides_guidel{guidel:02}.fa')
        logging.info(basename(guidesfap))
        if not exists(guidesfap):
            with open(guidesfap,'w') as f:
                for gi in dguides.index:
                    f.write('>{}\n{}\n'.format(gi.replace(' ','_'),dguides.loc[gi,'guide+PAM sequence']))

        ## BWA alignment command is adapted from cripror 
        ## https://github.com/rraadd88/crisporWebsite/blob/master/crispor.py
        # BWA: allow up to X mismatches
        # maximum number of occurences in the genome to get flagged as repeats. 
        # This is used in bwa samse, when converting the sam file
        # and for warnings in the table output.
        mismatches_max = 2

        MAXOCC = 60000

        # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
        MFAC = 2000000/MAXOCC
        genomed = dirname(genomep)

        bwaM = MFAC*MAXOCC # -m is queue size in bwa

        guidessap = os.path.join(datatmpd,f"01_guides_guidel{guidel:02}.sa" )

        logging.info(basename(guidessap))
        if not exists(guidessap):
            cmd=f"bwa aln -t 1 -o 0 -m {bwaM} -n {mismatches_max} -k {mismatches_max} -N -l {guidel} {genomep} {guidesfap} > {guidessap} 2> {guidessap}.log"
            print(cmd)
            os.system(cmd)

        guidessamp = os.path.join(datatmpd, f'01_guides_guidel{guidel:02}.sam')  
        logging.info(basename(guidessamp))        
        if not exists(guidessamp):
            cmd=f"bwa samse -n {MAXOCC} {genomep} {guidessap} {guidesfap} > {guidessamp} 2> {guidessamp}.log"
            print(cmd)
            os.system(cmd)

        guidessamp = os.path.join(datatmpd, f'01_guides_guidel{guidel:02}.sam')
        logging.info(basename(guidessamp))
        if not exists(guidessamp):
            cmd=f"bwa samse -n {MAXOCC} {genomep} {guidessap} {guidesfap} > {guidessamp} 2> {guidessamp}.log"
            print(cmd)
            os.system(cmd)  
            
def guidessam2dalignbed(datatmpd):
    """
    Processes SAM file to get the genomic coordinates in BED format
    step#2

    :param cfg: configuration dict
    """
    
    bed_colns = ['chromosome','start','end','id','NM','strand']

    
    alignmentbedp = os.path.join(datatmpd, '02_alignment.bed')
    dalignbedp = os.path.join(datatmpd, '02_dalignbed.tsv')
    logging.info(basename(dalignbedp))
    if not exists(alignmentbedp):
        guidessamps=glob(f'{datatmpd}/01_guides_guidel*.sam')
        for guidessamp in guidessamps:
            if os.stat(guidessamp).st_size != 0:
                samfile=pysam.AlignmentFile(guidessamp, "rb")
                dalignbed=pd.DataFrame(columns=bed_colns)
                for read in samfile.fetch():
                    algnids=[]
                    algnids.append('{}|{}{}|{}|{}'.format(read.reference_name,
                                 ('-' if read.is_reverse else '+'),read.positions[0],read.cigarstring,read.get_tag('NM')))
                    if read.has_tag('XA'):
                            algnids+=['|'.join(s.split(',')) for s in read.get_tag('XA').split(';') if len(s.split(','))>1]
                    chroms=[]
                    starts=[]
                    ends=[]
                    ids=algnids
                    NMs=[]
                    strands=[]
                    for a in algnids:
                        strand=a.split('|')[1][0]
                        chroms.append(a.split('|')[0])
                        if strand=='+':
                            offset=0
                        elif strand=='-':
                            offset=0                    
                        starts.append(int(a.split('|')[1][1:])+offset)
                        ends.append(int(a.split('|')[1][1:])+str2num(a.split('|')[2])+offset)
                        NMs.append(a.split('|')[3])
                        strands.append(strand)
                        del strand,offset
                    col2dalignbed={'chromosome':chroms,
                                       'start':starts,
                                       'end':ends,
                                       'id':ids,
                                       'NM':NMs,
                                       'strand':strands}
                    dalignbed_=pd.DataFrame(col2dalignbed)
                    dalignbed_['guide: id']=read.qname#.replace('_',' ')
                    dalignbed = dalignbed.append(dalignbed_,ignore_index=True,sort=True)
                samfile.close()
            else:
                logging.warning(f"file is empty: {guidessamp}")
            dalignbed.to_csv(dalignbedp,sep='\t')

            dalignbed['chromosome']=dalignbed.apply(lambda x : str2numorstr(x['chromosome']),axis=1)
            dalignbed=dalignbed.sort_values(['chromosome','start','end'], ascending=[True, True, True])
            dalignbed.loc[:,bed_colns].to_csv(alignmentbedp,sep='\t',
                            header=False,index=False,
                            chunksize=5000)
                        
def dalignbed2annotationsbed(datatmpd, alignmentbedp, genomegffp):
    """
    Get annotations from the aligned BED file
    step#3

    :param cfg: configuration dict
    """
    
    alignmentbedsortedp = alignmentbedp + '.sorted.bed'
    logging.info(basename(alignmentbedsortedp))
    if not exists(alignmentbedsortedp):
        cmd='{} sort -i {} > {}'.format('bedtools',alignmentbedp,alignmentbedsortedp)
        os.system(cmd)
    
    genomegffsortedp = genomegffp + '.sorted.gff3.gz'
    logging.info(basename(genomegffsortedp))
    if not exists(genomegffsortedp):    
        cmd=f"bedtools sort -i {genomegffp} > {genomegffsortedp}"
        os.system(cmd)

    annotationsbedp='{}/03_annotations.bed'.format(datatmpd)
    logging.info(basename(annotationsbedp))
    if not exists(annotationsbedp):    
        cmd=f"bedtools intersect -wa -wb -loj -a {alignmentbedsortedp} -b {genomegffsortedp} > {annotationsbedp}"
        os.system(cmd)
                
def dalignbed2dalignbedguides(datatmpd, dalignbedp, dguidesp):
    dalignbed=del_Unnamed(pd.read_csv(dalignbedp,sep='\t',keep_default_na=False))     #bed

    dguides=set_index(del_Unnamed(pd.read_csv(dguidesp,sep='\t',keep_default_na=False)),'guide: id') #pam tsv

    dalignbedguidesp = os.path.join(datatmpd ,"04_dalignbedguides.tsv")
    logging.info(basename(dalignbedguidesp))
    if not exists(dalignbedguidesp):
        dalignbed=pd.merge(dalignbed,dguides,on='guide: id',suffixes=('', '.1'))
        dalignbed.to_csv(dalignbedguidesp,'\t')
                
def alignmentbed2dalignedfasta(genomep ,datatmpd, alignmentbedp, dalignedfastap):

    logging.info(basename(dalignedfastap))
    if not exists(dalignedfastap):
        alignedfastap='{}/05_alignment.fa'.format(datatmpd)
        if not exists(alignedfastap):
            cmd=f"bedtools getfasta -s -name -fi {genomep} -bed {alignmentbedp} -fo {alignedfastap}"
            os.system(cmd)

        dalignedfasta=fa2df(alignedfastap)
        dalignedfasta.columns=['aligned sequence']
        dalignedfasta=dalignedfasta.loc[(dalignedfasta.apply(lambda x: not 'N' in x['aligned sequence'],axis=1)),:] #FIXME bwa aligns to NNNNNs
        dalignedfasta.index=[i.split('(')[0] for i in dalignedfasta.index] # for bedtools 2.27, the fasta header now has hanging (+) or (-)
        dalignedfasta.index.name='id'
        dalignedfasta.to_csv(dalignedfastap,sep='\t')
        
def dalignbed2dalignbedguidesseq(dalignbedguidesp, dalignedfastap, dalignbedguidesseqp):
    
    dalignbedguides = del_Unnamed(pd.read_csv(dalignbedguidesp, sep='\t', keep_default_na=False))
    dalignedfasta = del_Unnamed(pd.read_csv(dalignedfastap, sep='\t', keep_default_na=False))
    
    logging.info(basename(dalignbedguidesseqp))
    
    if not exists(dalignbedguidesseqp):        
        dalignbedguidesseq=pd.merge(dalignbedguides, dalignedfasta, on='id', suffixes=('', '.2'))
        dalignbedguidesseq=dalignbedguidesseq.dropna(subset=['aligned sequence'],axis=0)

        # dalignbed.index.name='id'
        dalignbedguidesseq=dalignbedguidesseq.drop_duplicates()
        dalignbedguidesseq.to_csv(dalignbedguidesseqp,sep='\t')

def dalignbedguidesseq2dalignbedstats(dalignbedguidesseqp, dalignbedstatsp):
    dalignbedguidesseq=del_Unnamed(pd.read_csv(dalignbedguidesseqp, sep='\t', keep_default_na=False))

    logging.info(basename(dalignbedstatsp))
    if not exists(dalignbedstatsp):
        df=dalignbedguidesseq.apply(lambda x: align(x['guide+PAM sequence'],x['aligned sequence']),
                               axis=1).apply(pd.Series)
        df.columns=['alignment','alignment: score']
        dalignbedstats=dalignbedguidesseq.join(df)
        del df
        dalignbedstats.to_csv(dalignbedstatsp,sep='\t')
                
def dannots2dalignbed2dannotsagg(datatmpd ,dannotsaggp, annotationsbedp):
    
    daannotp=f'{datatmpd}/08_dannot.tsv'

    logging.info(basename(daannotp))
    
    gff_colns = ['chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    bed_colns = ['chromosome','start','end','id','NM','strand']


    if ((not exists(daannotp)) and (not exists(dannotsaggp))):
            dannots=pd.read_csv(annotationsbedp,
                               sep='\t',
                               keep_default_na=False,
                               names= bed_colns + [c + ' annotation' if c in set(bed_colns).intersection(gff_colns) else c for c in gff_colns ],
                               low_memory=False)
            dannots=del_Unnamed(dannots)
            dannots=dannots.set_index('id')
            dannots['annotations count']=1

            # separate ids from attribute columns
            dannots=lambda2cols(dannots,lambdaf=gffatributes2ids,
                                in_coln='attributes',
                                to_colns=['gene name','gene id','transcript id','protein id','exon id'])
            dannots['annotation coordinate']=dannots.apply(lambda x: '{}:{}-{}({})'.format(x['chromosome annotation'],x['start annotation'],x['end annotation'], x['strand annotation']),axis=1)
            logging.debug('or this step takes more time?')
            dannots.to_csv(daannotp,sep='\t')
    else:
            dannots=pd.read_csv(daannotp,sep='\t',low_memory=False,keep_default_na=False)
            dannots=del_Unnamed(dannots)
    logging.info(basename(dannotsaggp))

    if not exists(dannotsaggp):
            if not 'dannots' in locals():
                dannots=pd.read_table(daannotp,low_memory=False,keep_default_na=False)
            dannots=del_Unnamed(dannots)
            dannots=dannots.reset_index()

            dannotsagg=pd.DataFrame(dannots.groupby('id')['annotations count'].agg('sum'))-1
            dannotsagg.loc[dannotsagg['annotations count']==0,'region']='intergenic'
            dannotsagg.loc[dannotsagg['annotations count']!=0,'region']='genic'

            alignids=dannots['id'].unique()#[:15]
            logging.debug('start of the slowest step')

            for alignidi in range(len(alignids)):
                alignid=alignids[alignidi]
                dannoti=dannots.loc[dannots['id']==alignid,:]
                if len(dannoti.shape)==1:
                    dannoti=pd.DataFrame(dannoti).T
                dannoti=dannoti.dropna(subset=['type'])
                if len(dannoti)!=0:
                    dannoti=dannoti.loc[dannoti['type']!='chromosome',:].drop_duplicates(subset=['start annotation','end annotation'])
                    for col in ['type','gene name','gene id','transcript id','protein id','exon id']:    
                        dannotsagg.loc[alignid,col+'s']=";".join(np.unique(dannoti[col].fillna('nan').tolist()))
            logging.debug('end of the slowest step')

            del dannots    
            dannotsagg=dannotsagg.reset_index()
            dannotsagg.to_csv(dannotsaggp,sep='\t')
            
def rescale(a,mn=None):
    a=(a-a.min())/(a.max()-a.min())
    if not mn is None:
        a=1-a
        a=a*(1-mn)
        a=1-a        
    return a

# wt = 'ATTCAATCCTATGCTGACTTGGG'
# off = 'ATTCAATCCTATGATGACTTGGG'
def get_cfdscore(wt,off):
    m_wt = re.search('[^ATCG]',wt)
    m_off = re.search('[^ATCG]',off)    
    if (m_wt is None) and (m_off is None):
        if len(wt)== 23 and len(off)== 23:
            pam = off[-2:]
            sg = off[:-3]
            cfd_score = calc_cfd(wt,sg,pam)
    #         print("CFD score: "+str(cfd_score))
            return cfd_score
    else:
        np.nan

def get_ppamdist(guidelength,pamlength,pam_position,ppamdist_min,pmutatpam):
    """
    Get penalties set based on distances of the mismatch/es from PAM 

    :param guidelength: length of guide sequence
    :param pamlength: length of PAM sequence
    :param pam_position: PAM location 3' or 5'
    :param ppamdist_min: minimum penalty
    :param pmutatpam: penalty for mismatch at PAM
    """
    x=np.arange(pamlength)
#     ppamdist_pam=np.zeros(pamlength)
    ppamdist_pam=np.full([1, pamlength], pmutatpam)[0]
    
    guide_positions=np.arange(guidelength)
    guide_positions=rescale(guide_positions)    
    coeffs=np.array([ 0.86947513, -1.15561009, -0.0826878 ,  0.77644198])
    ppamdist=np.polyval(coeffs,guide_positions)
    ppamdist=rescale(ppamdist,mn=ppamdist_min)
    ppamdist=np.concatenate((ppamdist,ppamdist_pam))
#     ppamdist=1-ppamdist
    if pam_position=='up':
        ppamdist=ppamdist[::-1]
    elif pam_position=='down':
        ppamdist=ppamdist
    else:
        raise(ValueError('ppamdist is invalid: {ppamdist}'))
    return ppamdist

def get_beditorscore_per_alignment(NM,genic,alignment,pam_length,pam_position,
                    pentalty_genic=0.5,
                    pentalty_intergenic=0.9,
                    pentalty_dist_from_pam=0.1,
                    test=False,debug=False):
    """
    Calculates beditor score per alignment between guide and genomic DNA.

    :param NM: Hamming distance
    :param mismatches_max: Maximum mismatches allowed in alignment
    :param genic: True if guide aligns to genic regions, else (intergenic) False.
    :param alignment: Symbol '|' means a match, '.' means mismatch and ' ' means gap. e.g. |||||.||||||||||.||||.| 
    :param pentalty_genic: penalty for genic alignment
    :param pentalty_intergenic: penalty for intergenic alignment
    :param pentalty_dist_from_pam: maximum pentalty for a mismatch at PAM () 
    :returns: beditor score per alignment.
    """
    pmutatpam=1e-300
    if not pd.isnull(alignment):
        if NM!=0:            
            pentalty_region_of_alignment=pentalty_genic if genic else pentalty_intergenic 

            mutations_penalties=np.array([0 if s=='|' else 1 for s in alignment])
            dist_from_pam_penalties=get_ppamdist(guidelength=len(alignment)-pam_length,
                                                 pamlength=pam_length,
                                                 pam_position=pam_position,
                                                 ppamdist_min=pentalty_dist_from_pam,
                                                pmutatpam=pmutatpam)
            mutations_penalties_multi=mutations_penalties*dist_from_pam_penalties
            mutations_penalties_multi=mutations_penalties_multi[mutations_penalties_multi != 0]
            if len(mutations_penalties_multi)==0:
                # all the penalties are 0
                penality_cum_dist_from_pam=0
            else:
                if pmutatpam in mutations_penalties_multi:
                    # mutation at pam
                    penality_cum_dist_from_pam=0
                else:
                    penality_cum_dist_from_pam=np.prod(mutations_penalties_multi)
            if test:
                print(dist_from_pam_penalties)
                print(mutations_penalties)
                print(mutations_penalties_multi)            
                print(penality_cum_dist_from_pam)            
            if debug:
                return pentalty_region_of_alignment,penality_cum_dist_from_pam                
            return pentalty_region_of_alignment*penality_cum_dist_from_pam
        else:
            if debug:
                return np.nan,np.nan                
            return 1
    else:
        if debug:
            return np.nan,np.nan                
        return np.nan

def dannotsagg2dannots2dalignbedannot(dannotsaggp, dalignbedstatsp, dalignbedannotp):
    dannotsagg=del_Unnamed(pd.read_csv(dannotsaggp, sep='\t', keep_default_na=False))
    dalignbedstats=del_Unnamed(pd.read_csv(dalignbedstatsp, sep='\t', keep_default_na=False))
    
    dalignbedannotp
    logging.info(basename(dalignbedannotp))
    
    if not exists(dalignbedannotp):
        
        dalignbedannot=dalignbedstats.set_index('id').join(set_index(dannotsagg,'id'),
                                              rsuffix=' annotation')
        dalignbedannot['NM']=dalignbedannot['NM'].apply(int)
        
        dalignbedannot['beditor score']=dalignbedannot.apply(lambda x : get_beditorscore_per_alignment(NM=x['NM'],
                               genic=True if x['region']=='genic' else False,
                               alignment=x['alignment'],
                               pam_length=len(x['PAM']),
                               pam_position=x['original position']
                                ),axis=1) 
        dalignbedannot['CFD score']=dalignbedannot.apply(lambda x : get_cfdscore(x['guide+PAM sequence'].upper(), 
                                                                                 x['aligned sequence'].upper()), axis=1)
        dalignbedannot['CFD score']=dalignbedannot['CFD score'].fillna(0)
        dalignbedannot.to_csv(dalignbedannotp,sep='\t')
                
def dalignbedannot2daggbyguide(datatmpd, dalignbedannotp, dbepamsp ,be_names, dofftargetsp, daggbyguidep):
    
        # dbepams
        cols_dbes=['distance of codon start from PAM: maximum',
         'distance of codon start from PAM: minimum',
         'distance of mutation from PAM: maximum',
         'distance of mutation from PAM: minimum',
         'method',
         'nucleotide',
         'nucleotide mutation',
         'strand']

        dalignbedannot=del_Unnamed(pd.read_csv(dalignbedannotp,sep='\t',low_memory=False,keep_default_na=False))

        # daggbyguidep=f'{datatmpd}/10_daggbyguide.tsv'      
        logging.info(basename(daggbyguidep))  

        if not exists(daggbyguidep):
            dbepams=pd.read_table(dbepamsp, keep_default_na=False)

            dBEs=dbepams.loc[:,cols_dbes]
            dBEs=dBEs.loc[dBEs['method'].isin(be_names),:]        
            daggbyguide=dalignbedannot.loc[(dalignbedannot['NM']==0),['guide: id','guide+PAM sequence','gene names', 'gene ids','transcript ids']].drop_duplicates(subset=['guide: id'])
            if len(daggbyguide)!=0:
                daggbyguide=set_index(daggbyguide,'guide: id')            
                guideids=daggbyguide.index.tolist()
                for gi in range(len(guideids)):
                    gid=guideids[gi]
                    dalignbedannoti=dalignbedannot.loc[dalignbedannot['guide: id']==gid,:]
                    if len(dalignbedannoti.shape)==1:
                        dalignbedannoti=pd.DataFrame(dalignbedannoti).T
                    for col in ['types','gene names','gene ids','transcript ids','protein ids','exon ids']:
                        if (col in daggbyguide) and (col in dalignbedannoti):
                            daggbyguide.loc[gid,col]=";".join(np.unique(dalignbedannoti[col].fillna('nan').tolist()))

                for guideid in daggbyguide.index:
                    dalignbedannotguide=dalignbedannot.loc[(dalignbedannot['guide: id']==guideid),:]
                    daggbyguide.loc[guideid,'beditor score']=get_beditorscore_per_guide(guide_seq=dalignbedannotguide['guide+PAM sequence'].unique()[0], 
                strategy=dalignbedannotguide['strategy'].unique()[0],
                align_seqs_scores=dalignbedannotguide['beditor score'],
                dBEs=dBEs
    #                                        test=cfg['test']
                )

                    daggbyguide.loc[guideid,'CFD score']=dalignbedannotguide['CFD score'].mean() #FIXME if mean is not appropriate
                daggbyguide['beditor score (log10)']=daggbyguide['beditor score'].apply(np.log10)
                dalignbedannot['alternate alignments count']=1
                daggbyguide=daggbyguide.join(pd.DataFrame(dalignbedannot.groupby('guide: id')['alternate alignments count'].agg('sum')))
                daggbyguide.to_csv(daggbyguidep,sep='\t')
                daggbyguide.to_csv(dofftargetsp,sep='\t')

def get_dguides2offtargets(genomep, genomegffp, datatmpd, dguidesp, dofftargetsp  ,dbepamsp):

    dguides= del_Unnamed(pd.read_csv(dguidesp,sep='\t',keep_default_na=False))
    if len(dguides)==0:
        logging.warning(f"dguides is empty.")

    step2doutp={
        1:'01_guides_guidel*.fa',
        2:'02_dalignbed.tsv',
        3:'03_annotations.bed',
        4:'04_dalignbedguides.tsv',
        5:'05_dalignedfasta.tsv',
        6:'06_dalignbedguidesseq.tsv',
        7:'07_dalignbedstats.tsv',
        8:'08_dannotsagg.tsv',
        9:'09_dalignbedannot.tsv',
        10:'10_daggbyguide.tsv',
        }

    dguidesp = dguidesp
    alignmentbedp = os.path.join(datatmpd, '02_alignment.bed')  
    dalignbedp = os.path.join(datatmpd, '02_dalignbed.tsv')    
    dalignbedguidesp = os.path.join(datatmpd, '04_dalignbedguides.tsv')  
    dalignedfastap = os.path.join(datatmpd, '05_dalignedfasta.tsv')  
    dalignbedguidesseqp = os.path.join(datatmpd, '06_dalignbedguidesseq.tsv')  
    dalignbedstatsp = os.path.join(datatmpd, '07_dalignbedstats.tsv') 
    dannotsaggp = os.path.join(datatmpd, '08_dannotsagg.tsv')  
    dalignbedannotp = os.path.join(datatmpd, '09_dalignbedannot.tsv') 
    daggbyguidep = os.path.join(datatmpd, '10_daggbyguide.tsv')  


    # dofftargetsp = os.path.join(basedir,'dofftargets.tsv')





    if not exists(dofftargetsp):
        # genomep = '/home/yanghe/githubcode/beditor/super_beditor/pub/GCF_11111/fasta/yb01/dna/genome.fna'
        dguides2guidessam(datatmpd,genomep,dguides)
        
        guidessam2dalignbed(datatmpd)
        
        # genomegffp = '/home/yanghe/githubcode/beditor/super_beditor/pub/GCF_11111/gff3/yb01/genome.gff'
        dalignbed2annotationsbed(datatmpd, alignmentbedp, genomegffp)   
        
        dalignbed2dalignbedguides(datatmpd, dalignbedp, dguidesp)  
        
        alignmentbed2dalignedfasta(genomep, datatmpd, alignmentbedp, dalignedfastap)
        
        dalignbed2dalignbedguidesseq(dalignbedguidesp, dalignedfastap, dalignbedguidesseqp)
        
        dalignbedguidesseq2dalignbedstats(dalignbedguidesseqp, dalignbedstatsp)

        annotationsbedp='{}/03_annotations.bed'.format(datatmpd)
        dannots2dalignbed2dannotsagg(datatmpd ,dannotsaggp, annotationsbedp)
        
        dannotsagg2dannots2dalignbedannot(dannotsaggp, dalignbedstatsp, dalignbedannotp)
        
        # dbepamsp = '/home/yanghe/githubcode/beditor/sgRNA_primer/input_output/YB01/00_input/dbepams.tsv'
        be_names = ['Target-AID']
        dalignbedannot2daggbyguide(datatmpd, dalignbedannotp, dbepamsp, be_names, dofftargetsp, daggbyguidep)
        
        import gc
        gc.collect() 






def main():

    p_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))

    basedir = os.path.join(p_dir,'data/YB01/04_offtargets/')
    if not exists(basedir):
        os.makedirs(basedir)

    dofftargetsp = os.path.join(basedir, 'dofftargets.tsv')
    dguidesp = os.path.join(p_dir, 'data/YB01/03_guides/','dguides.tsv')

    datatmpd = os.path.join(basedir, 'tmp/')

    if not exists(datatmpd):
        os.makedirs(datatmpd)

    genomep =  '/home/yanghe/baseEditor/beditor/pub/GCF001/fasta/yb01/dna/genome.fna'
    genomegffp =  '/home/yanghe/baseEditor/beditor/pub/GCF001/gff3/yb01/genome.gff'
    dbepamsp =  '/home/yanghe/baseEditor/beditor/data/yb01/00_input/dbepams.tsv'

    get_dguides2offtargets(genomep, genomegffp, datatmpd, dguidesp, dofftargetsp, dbepamsp)
    

   


if __name__ == '__main__':

    main()