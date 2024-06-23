import os
from os.path import exists,abspath,dirname
import pandas as pd
import logging
import numpy as np
import json
from os import makedirs


cols_dpam=['PAM', 'PAM position', 'guide length']


def reverse_complement_multintseq(seq,nt2complement):
    complement=[]
    for s in list(seq):
        for ss in nt2complement:
            if ss==s:
#                 print(nt2complement[s],s)
                complement.append(nt2complement[s])
                break
    return "".join(complement[::-1]    )


def get_nt2complement(): 
    nt2complement={'A':'T',
                  'G':'C',
                  'N':'N',
                  'R':'Y',
                  'S':'W',
                  'K':'M',
                   'B':'b',
                   'D':'d',
                   'H':'h',
                   'V':'v',
                   'N':'N',
                  }
    nt2complement.update(dict(zip(nt2complement.values(),nt2complement.keys())))
    return nt2complement
nt2complement=get_nt2complement()


multint2reg={'R':'[AG]',
'Y':'[CT]',
'S':'[GC]',
'W':'[AT]',
'K':'[GT]',
'M':'[AC]',
'B':'[CGT]',
'D':'[AGT]',
'H':'[ACT]',
'V':'[ACG]',
'N':'[ATGC]',}
multint2regcomplement={'R':'[TC]',
'Y':'[GA]',
'S':'[GC]',
'W':'[AT]',
'K':'[CA]',
'M':'[TG]',
'B':'[^A]',
'D':'[^C]',
'H':'[^G]',
'V':'[^T]',
'N':'[ATGC]',}

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

def df2info(df,col_searches=None):
    if len(df.columns)>5:
        print('**COLS**: ',df.columns.tolist())
    print('**HEAD**: ',df.loc[:,df.columns[:5]].head())
    print('**SHAPE**: ',df.shape)
    if not col_searches is None:
        cols_searched=[c2 for c1 in col_searches for c2 in df if c1 in c2]
        print('**SEARCHEDCOLS**:\n',cols_searched)
        print('**HEAD**: ',df.loc[:,cols_searched].head()) 


def del_Unnamed(df):
    """
    Deletes all the unnamed columns

    :param df: pandas dataframe
    """
    cols_del=[c for c in df.columns if 'Unnamed' in c]
    return df.drop(cols_del,axis=1)

def to_table_pqt(df,p):
    if not exists(dirname(p)) and dirname(p)!='':
        makedirs(dirname(p),exist_ok=True)
    df.to_parquet(p,engine='fastparquet',compression='gzip',)

def read_table_pqt(p):
    return del_Unnamed(pd.read_parquet(p,engine='fastparquet'))

def s2re(s,ss2re):
    for ss in ss2re:
        s=s.replace(ss,ss2re[ss])
    return s

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
        
def to_table(df,p):
    if p.endswith('.tsv') or p.endswith('.tab'):
        if not exists(dirname(p)) and dirname(p)!='':
            makedirs(dirname(p),exist_ok=True)
        df.to_csv(p,sep='\t')
    elif p.endswith('.pqt') or p.endswith('.parquet'):
        to_table_pqt(df,p)
    else: 
        logging.error(f'unknown extension {p}')      
        
def saveemptytable(step,doutp=None):
   
    dout=pd.DataFrame(columns=stepi2cols[step])
    logging.warning(f"saved enmpty table {doutp}")
    if not doutp is None:
        dout.to_csv(doutp,sep='\t')
    else:
        return dout 


def reverse_complement_multintseqreg(seq,multint2regcomplement,nt2complement):    
    complement=[]  
    for s in list(seq):  
        if s in multint2regcomplement.keys():
            for ss in multint2regcomplement:
                if ss==s:   
    #                 print(nt2complement[s],s)
                    complement.append(multint2regcomplement[s])
                    break
        elif s in nt2complement.keys():
            for ss in nt2complement:
                if ss==s:
                    complement.append(nt2complement[s])
                    break            
        else:
            logging.error(f'odd character {s} in seq {seq}')
        
    return "".join(complement[::-1] )  


def dinnucleotide2dsequencesproper(dsequences,dmutagenesis,dbug=False):
    """
    Makes dseqeunces dataframe of nucleotide mutation format compatible to guide design modules   使核苷酸突变格式的dseqeunces数据框架与指导设计模块兼容

    :param dsequences: dsequences dataframe
    :param dmutagenesis: dmutagenesis dataframe
    """

    dmutagenesis=dmutagenesis.loc[(dmutagenesis['position of mutation in codon']==2),:]
    # print(dsequences['nucleotide wild-type','nucleotide mutation'])   

    dsequences=pd.merge(dsequences,dmutagenesis,      
             left_on=['nucleotide wild-type','nucleotide mutation','codon: wild-type', 'codon: mutation'],
             right_on=['nucleotide: wild-type','nucleotide: mutation','codon','codon mutation'],how='inner',
             suffixes=['', ': dmutagenesis'])  
  

    if len(dsequences)!=0:
        dsequences['transcript: id']=dsequences['genome coordinate']
        dsequences['aminoacid mutation']=dsequences['amino acid mutation']
        dsequences['aminoacid: wild-type']=dsequences['amino acid']
        if dbug:
            dsequences['codon: wild-type']=dsequences['codon']
            df2info(dsequences,'nucle')
            df2info(dmutagenesis,'nucle')
        dsequences['id']=dsequences.apply(lambda x: f"{x['genome coordinate']}|{x['method']}|{x['mutation on strand'].replace(' strand','')}|{x['nucleotide wild-type']}:{x['nucleotide mutation']}|{x['codon: wild-type']}:{x['codon mutation']}",axis=1)
    else:
        logging.warning('empty dsequences after merging with dmutagenesis')
    cols_missing=[c for c in stepi2cols[1] if not c in dsequences]
    for c in cols_missing:
        dsequences[c]=np.nan
    return dsequences,dmutagenesis


def dpam2dpam_strands(dpam,pams):
    """
    Duplicates dpam dataframe to be compatible for searching PAMs on - strand

    :param dpam: dataframe with pam information
    :param pams: pams to be used for actual designing of guides.
    """
    
    dpam=del_Unnamed(dpam)
    dpam['rPAM']=dpam.apply(lambda x : s2re(x['PAM'],multint2reg) ,axis=1)
    dpam=set_index(dpam,'PAM')
    dpam['strand']='+'
    dpamr=pd.DataFrame(columns=dpam.columns)
    dpam.loc[:,'reverse complement']=np.nan
    dpam.loc[:,'original']=np.nan
    for pam in dpam.index:
        pamr=reverse_complement_multintseq(pam,nt2complement)
        dpam.loc[pam,'reverse complement']=pamr
        dpam.loc[pam,'original']=pam
        dpamr.loc[pamr,'original']=pam 
        dpam.loc[pam,'original position']=dpam.loc[pam,'PAM position']
        dpamr.loc[pamr,'original position']=dpam.loc[pam,'PAM position']
        dpamr.loc[pamr,['PAM position','guide length']]=dpam.loc[pam,['PAM position','guide length']]
        dpamr.loc[pamr,['rPAM']]=reverse_complement_multintseqreg(pam,multint2regcomplement,nt2complement)    
    dpamr['PAM position']= dpamr.apply(lambda x: 'up' if x['PAM position']=='down' else 'down',axis=1)
    dpamr['strand']='-'
    dpam_strands=dpam.append(dpamr,sort=True)
    dpam_strands.index.name='PAM'
    dpam_strands.loc[:,'is a reverse complement']=pd.isnull(dpam_strands.loc[:,'reverse complement'])
    pams_strands=pams+dpam_strands.loc[pams,'reverse complement'].dropna().tolist()
    dpam_strands=dpam_strands.loc[pams_strands,:]
    return dpam_strands


def get_be2dpam(din):
    """
    make BE to dpam mapping i.e. dict
    
    :param din: df with BE and PAM info all cols_dpam needed
    """
    be2dpam={}
    be2pam=din.loc[:,['method','PAM']].drop_duplicates().set_index('method').to_dict()['PAM']
#     if test:
#         print(be2pam)
    for be in be2pam:
        pam=be2pam[be]
        dpam=din.loc[((din['PAM']==pam) & (din['method']==be) & (din['strand']=='+')),cols_dpam]
        dpam_strands=dpam2dpam_strands(dpam,pams=[pam])
        be2dpam[be]=set_index(dpam_strands,'PAM')                
    return be2dpam


def get_pam_searches(dpam,seq,pos_codon):
    """
    Search PAM occurance

    :param dpam: dataframe with PAM sequences
    :param seq: target sequence
    :param pos_codon: reading frame
    :param test: debug mode on
    :returns dpam_searches: dataframe with positions of pams
    """
    import regex as re
    def get_guide_pam(match,pam_stream,guidel,pos_codon):  
        if pam_stream=='down':
            # IIIIIIINGG
            # 0123456789
            seq_guidepam=seq[match.span()[0]-guidel:match.span()[1]]
            seq_guide=seq[match.span()[0]-guidel:match.span()[0]]
            seq_pam=seq[match.span()[0]:match.span()[1]]   
            dist_codon=pos_codon-match.span()[0]
#             if test:
#                 print(match.span()[0]-pos_codon)
        elif pam_stream=='up':
            # TTNIIIIIII
            # 0123456789
            seq_guidepam=seq[match.span()[0]:match.span()[1]+guidel]
            seq_guide=seq[match.span()[1]:match.span()[1]+guidel]
            seq_pam=seq[match.span()[0]:match.span()[1]]
            dist_codon=pos_codon-match.span()[1]
        if seq_pam!=match.group():
            logging.error(f'indexing is wrong:seq_guidepam: {seq_guidepam}, seq_guide: {seq_guide}, seq_pam: {seq_pam},match.group(): {match.group()}')               
        return seq_guidepam,seq_guide,seq_pam,abs(dist_codon)
    pams=dpam.index.tolist()
    dpamposs=pd.DataFrame(columns=['guide+PAM sequence','guide sequence','PAM sequence'])
    pamposi=0
#     print(dpam)
    for pam in pams:
#         print(pam,pams)
        matchesiter=re.finditer(dpam.loc[pam,'rPAM'], seq, overlapped=True)
        for match in matchesiter:
            dpamposs.loc[pamposi,'position of PAM ini'],dpamposs.loc[pamposi,'position of PAM end']=match.span()   #返回一个元组
            dpamposs.loc[pamposi,'position of PAM end']=dpamposs.loc[pamposi,'position of PAM end']-1                
            dpamposs.loc[pamposi,'guide+PAM sequence'],dpamposs.loc[pamposi,'guide sequence'],dpamposs.loc[pamposi,'PAM sequence'],dpamposs.loc[pamposi,'distance of codon from PAM']\
            =get_guide_pam(match,dpam.loc[pam,'PAM position'],dpam.loc[pam,'guide length'],pos_codon)
            dpamposs.loc[pamposi,'PAM']=pam
            pamposi+=1
    dpamposs['codon: from pam search']=seq[pos_codon:pos_codon+3]
    dpamposs['guide sequence']=dpamposs['guide sequence'].fillna('')
    if len(dpamposs)==0:
        return None
    dpamposs['guide sequence length']=dpamposs.apply(lambda x : len(x['guide sequence']),axis=1)

    dpamposs['franking sequence']=seq
    dpamposs=dpam.join(set_index(dpamposs,'PAM'),how='right')
    dpamposs.index.name='PAM'
    return dpamposs


def guide2dpositions(x,dbug=False): 
    """
    Get positions of guides relative to the target site and PAM sequence
    Note:
    Index and flank sequence based indexing are 0-based
    Distances and positions from pam are 1-based

    :param x: lambda section of dguides dataframe  
    """
    dpositions=pd.DataFrame(index=range(45),
                           columns=['guide+PAM sequence'])
    dpositions.index.name='PAM position'

    dpositions.loc[x['position of PAM ini']:x['position of PAM end'],'location PAM']=True
    if x['PAM position']=='up':
        dpositions.loc[x['position of PAM end']+1:x['position of PAM end']+x['guide sequence length'],'location guide']=True
        dpositions.loc[x['position of PAM end']+x['distance of mutation from PAM: minimum']:x['position of PAM end']+x['distance of mutation from PAM: maximum'],'location window']=True        
    elif x['PAM position']=='down':
        dpositions.loc[x['position of PAM ini']-x['guide sequence length']:x['position of PAM ini']-1,'location guide']=True
        dpositions.loc[x['position of PAM ini']-x['distance of mutation from PAM: maximum']:x['position of PAM ini']-x['distance of mutation from PAM: minimum'],'location window']=True        

    dpositions.loc[21:23,'location codon']=True
    dpositions.loc[21-1+x['position of mutation in codon'],'location mutation']=True
    dpositions[[c for c in dpositions if 'location' in c]]=dpositions[[c for c in dpositions if 'location' in c]].fillna(False)
    dpositions['location guide+PAM']=dpositions['location guide'] | dpositions['location PAM']
    if x['original position']=='down' and x['strand']=='+':
        dpositions['position from PAM']=np.array(dpositions.index.tolist())-dpositions.loc[dpositions['location PAM'],:].index.tolist()[0]
    elif x['original position']=='down' and  x['strand']=='-':
        dpositions['position from PAM']=np.array(dpositions.index.tolist())-dpositions.loc[dpositions['location PAM'],:].index.tolist()[-1]
        dpositions['position from PAM']=dpositions['position from PAM']*-1    
    elif x['original position']=='up' and x['strand']=='+':
        dpositions['position from PAM']=dpositions.loc[dpositions['location PAM'],:].index.tolist()[-1]-np.array(dpositions.index.tolist())
        dpositions['position from PAM']=dpositions['position from PAM']*-1    
    elif x['original position']=='up' and  x['strand']=='-':
        dpositions['position from PAM']=np.array(dpositions.index.tolist())-dpositions.loc[dpositions['location PAM'],:].index.tolist()[0]
        dpositions['position from PAM']=dpositions['position from PAM']*-1    

    dpositions.loc[dpositions['location mutation'],'nucleotide wild-type']=x['nucleotide']
    dpositions.loc[dpositions['location mutation'],'nucleotide mutation']=x['nucleotide mutation']

    dpositions.loc[dpositions['location codon'],'codon wild-type']=list(x['codon: wild-type'])
    dpositions.loc[dpositions['location codon'],'codon mutation']=list(x['codon mutation'])
    if x['strand']=='+':
        dpositions.loc[dpositions['location guide+PAM'],'guide+PAM sequence']=list(x['guide+PAM sequence'])
        activity_sequence=''.join(dpositions.loc[dpositions['location window'],'guide+PAM sequence'].tolist())
    elif x['strand']=='-':
        dpositions.loc[dpositions['location guide+PAM'],'guide+PAM sequence']=list(x['guide+PAM sequence'])[::-1]
        activity_sequence=''.join(dpositions.loc[dpositions['location window'],'guide+PAM sequence'].tolist())[::-1]
    posmut=dpositions.loc[dpositions['location mutation'],:].index[0]
    posmutfrompam=int(dpositions.loc[dpositions['location mutation'],'position from PAM'])
    distmutfrompam=abs(posmutfrompam)
    posguideini=dpositions.loc[dpositions['location guide'],:].index.min()
    posguideend=dpositions.loc[dpositions['location guide'],:].index.max()
    if dbug:
        print(x[['strategy','strand','distance of mutation from PAM']+[s for s in x.index if ('sequence' in s) or ('length' in s)]])
        print({'posmut':posmut,
               'posmutfrompam':posmutfrompam,
               'distmutfrompam':distmutfrompam,
               'posguideini':posguideini,
               'posguideend':posguideend,
              'activity_sequence':activity_sequence})
        return dpositions
    else:
        return posmut,posmutfrompam,distmutfrompam,posguideini,posguideend,activity_sequence
    

def make_guides(dseq, dmutagenesis, flankaas, dbepamsp):

    dseq=dseq.reset_index()
    dseq.index=range(len(dseq))
    if not 'pos control' in dseq:
        dseq['pos control']=False

    dseq_cols=['transcript: id','aminoacid: position','aminoacid: wild-type','codon: wild-type','id','pos control'] 
    dbepams= del_Unnamed(pd.read_table(dbepamsp ,keep_default_na=False) )
    be2dpam=get_be2dpam(dbepams)

    gierrfltmutpos=[]
    gierrdenan=[]
    gierrfltguidel=[]
    gierrpamnotfound=[]
    gierrcannotmutate=[]  

    for gi in dseq.index:
            for be in be2dpam:
                dmutagenesis_be=dmutagenesis.loc[dmutagenesis['method']==be,:]
                if len(dmutagenesis_be)==0:
                    continue
    #             if cfg['mutations']=='mutations':
                dseqi=pd.DataFrame(dseq.loc[gi,dseq_cols+['amino acid mutation']]).T
                dmutagenesis_gi=pd.merge(dseqi,
                    dmutagenesis_be,
                    how='inner',
                    left_on=['aminoacid: wild-type','codon: wild-type','amino acid mutation'],
                    right_on=['amino acid','codon','amino acid mutation'])     
                if len(dmutagenesis_gi)!=0:
                    pos_codon=(flankaas)*3
                    dpam=be2dpam[be]
                    dpamsearches=get_pam_searches(dpam=dpam,
                         seq=dseq.loc[gi,'transcript: sequence'],
                         pos_codon=pos_codon)
                    if dpamsearches is None:
                        continue
                    if len(dpamsearches)!=0:
                         # filter by guide length
                        dpamsearchesflt=dpamsearches.loc[dpamsearches['guide length']==dpamsearches['guide sequence length'],:]
                        dpamsearches_strategy=pd.merge(dpamsearchesflt.reset_index(),dmutagenesis_gi.reset_index(),
                                         how='inner',
                                         on=['strand'],suffixes=['',': dmutagenesis_gi'])
                        if len(dpamsearches_strategy)!=0:                                 
                            if not 'dguides' in locals():
                                dguides=dpamsearches_strategy.copy()
                            else:
                                dguides=dguides.append(dpamsearches_strategy)
                                del dpamsearches_strategy
                        else:
                            gierrdenan.append(gi)
                    else:
                        gierrpamnotfound.append(gi)
                else:
                    gierrcannotmutate.append(gi)

    gierrfltmutpos=[]
    gierrdenan=[]
    gierrfltguidel=[]
    gierrpamnotfound=[]
    gierrcannotmutate=[]

    err2idxs={'gierrfltmutpos':gierrfltmutpos,
                  'gierrdenan':gierrdenan,
                  'gierrfltguidel':gierrfltguidel,
                  'gierrpamnotfound':gierrpamnotfound,
                  'gierrcannotmutate':gierrcannotmutate,
                 }  




    if 'dguides' in locals():
        logging.info('#reverse complement guides on negative strand sequences')
        dguides.loc[:,'PAM']=dguides.apply(lambda x : reverse_complement_multintseq(x['PAM'],nt2complement) if x['is a reverse complement'] else x['PAM'],axis=1)
        for colseq in ['guide+PAM sequence','guide sequence','PAM sequence']:
            dguides.loc[:,colseq]=dguides.apply(lambda x : str(str2seq(x[colseq]).reverse_complement()) if x['is a reverse complement'] else x[colseq],axis=1)

        logging.info('get dposition')
        dpositions=dguides.apply(lambda x: guide2dpositions(x),axis=1).apply(pd.Series)
        dpositions.columns=['position of mutation','position of mutation from PAM',
                                'distance of mutation from PAM',
                               'position of guide ini','position of guide end','activity sequence']
        for col in dpositions:
            dguides[col]=dpositions[col]

        logging.info('filter by # of editable nts in activity seq')
        logging.info(dguides.shape)
        dguides_noflt=dguides.copy()        

        if len(dguides)!=0:
            dguides.loc[:,'strategy']=dguides.apply(lambda x: f"{x['method']};{x['strand']};@{int(x['distance of mutation from PAM'])};{x['PAM']};{x['codon']}:{x['codon mutation']};{x['amino acid']}:{x['amino acid mutation']};",axis=1)
            dguides.loc[:,'guide: id']=dguides.apply(lambda x: f"{x['id']}|{int(x['aminoacid: position']) if not pd.isnull(x['aminoacid: position']) else 'nucleotide'}|({x['strategy']})",axis=1)
            dguides.loc[:,'guide+PAM length']=dguides.apply(lambda x: len(x['guide+PAM sequence']),axis=1)
            dguides=dguides.drop_duplicates(subset=['guide: id'])

            logging.info('#filter by location of mutation within guide')
            dguides_neg_control=dguides.loc[dguides.apply(lambda x : False if (x['distance of mutation from PAM: minimum']<=abs(x['distance of mutation from PAM'])<=x['distance of mutation from PAM: maximum']) else True,axis=1),:]
            logging.info(dguides.shape)
            if len(dguides_neg_control)==0:
                dguides_neg_control=None
            dguides=dguides.loc[dguides.apply(lambda x : True if (x['distance of mutation from PAM: minimum']<=abs(x['distance of mutation from PAM'])<=x['distance of mutation from PAM: maximum']) else False,axis=1),:]            
            logging.info(f"filter by distance of mutation from PAM: dguides len :{len(dguides)}")   
            if len(dguides)!=0:
                logging.info(dguides.shape)
                dguides_pos_control=dguides.loc[dguides['pos control'],:]      
                if len(dguides_pos_control)==0:   
                    dguides_pos_control=None    
                dguides=dguides.loc[~dguides['pos control'],:]                
                logging.info(dguides['pos control'].sum())  
                logging.info(dguides.shape)


                return dguides,dguides_noflt,err2idxs,dguides_neg_control,dguides_pos_control
            else:
                return None,dguides_noflt,None,None,None       
        else:
            return None,dguides_noflt,None,None,None
    else:
        return None,None,None,None,None 
    

def get_guides( dguideslinp, dguides_nofltp, dmutagenesisp_03, dmutagenesisp_02,dsequences_01, dsequences_03, mutation_format, dbepamsp, flankaas, make_control_pos, make_control_neg ):

    if not exists(dguideslinp): 
        dmutagenesis= del_Unnamed(pd.read_csv(dmutagenesisp_02, sep='\t', keep_default_na=False))
        if mutation_format == 'nucleotide':
            dsequences = del_Unnamed(pd.read_csv(dsequences_01,sep='\t',keep_default_na=False))
            dsequences,dmutagenesis=dinnucleotide2dsequencesproper(dsequences,dmutagenesis)
        elif mutation_format == 'aminoacid':
            pass
        dsequences.to_csv(dsequences_03, sep='\t') 

        if not (len(dsequences)==0 or len(dmutagenesis)==0): 
            dmutagenesis['strand']=dmutagenesis.apply(lambda x : x['mutation on strand'].replace(' strand',''),axis=1) 
            dmutagenesis.to_csv(dmutagenesisp_03,sep='\t')
            dguideslin,dguides_noflt,err2idxs,dguides_neg_control,dguides_pos_control=make_guides(dsequences,dmutagenesis,flankaas,dbepamsp)
            if not dguides_noflt is None:
                dguides_noflt.to_csv(dguides_nofltp,sep='\t')
            if not ((dguideslin is None) and (err2idxs is None)):
                dguideslin.to_csv(dguideslinp,sep='\t')
                with open(dguideslinp+'.err.json', 'w') as f:
                    json.dump(err2idxs, f)
                if make_control_pos and not dguides_pos_control is None:
                    to_table(dguides_pos_control,f"{dguideslinp}.pos_control.tsv")
                if make_control_neg and not dguides_neg_control is None:
                    to_table(dguides_neg_control,f"{dguideslinp}.neg_control.tsv")
            else:
                logging.warning('no guides designed; saving an empty table.')
                saveemptytable(3,dguideslinp)   
        else:
            logging.warning('no guides designed; saving an empty table.')
            saveemptytable(3,dguideslinp)

    import gc
    gc.collect 




def main():

    p_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
    basedir = os.path.join(p_dir, 'data/YB01/03_guides/')


    if not exists(basedir):
        os.makedirs(basedir)
    dguideslinp = os.path.join(basedir, 'dguides.tsv')      
    dguides_nofltp = os.path.join(basedir, 'dguides_noflt.tsv')
    dmutagenesisp_03 = os.path.join(basedir, 'dmutagenesis.tsv')

    dmutagenesisp_02 = os.path.join(p_dir, 'data/YB01/02_mutagenesis/', 'dmutagenesis.tsv')
    dsequences_01 = os.path.join(p_dir, 'data/YB01/01_sequences/', 'dsequences.tsv')
    dsequences_03 = os.path.join(p_dir, 'data/YB01/03_guides/', 'dsequences.tsv')

    mutation_format = 'nucleotide'    
    dbepamsp = '/home/yanghe/githubcode/beditor/sgRNA_primer/input_output/YB01/00_input/dbepams.tsv'
    flankaas=7
    make_control_pos = None
    make_control_neg = None


    get_guides( dguideslinp, dguides_nofltp, dmutagenesisp_03, dmutagenesisp_02,dsequences_01, dsequences_03, mutation_format, dbepamsp, flankaas, make_control_pos, make_control_neg)
    
    

if __name__ == '__main__':

    main()  