import os
from os.path import exists,abspath,dirname,realpath
import pandas as pd
import logging
from Bio import SeqIO, Alphabet, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO
from Bio import Data

aminoacids=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","*"]
cols_dbes=['distance of codon start from PAM: maximum',
    'distance of codon start from PAM: minimum',
    'distance of mutation from PAM: maximum',
    'distance of mutation from PAM: minimum',
    'method',
    'nucleotide',
    'nucleotide mutation',
    'strand']


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

cuspp = os.path.join ( dirname( dirname( realpath(__file__) ) ), "input/data/64_1_1_all_nuclear.cusp.txt") 

# cuspp = "/hpcfs/fhome/yangchh/editSeqDesign/self/AutoBeditor/input/data/64_1_1_all_nuclear.cusp.txt"

def del_Unnamed(df):
    """
    Deletes all the unnamed columns

    :param df: pandas dataframe
    """
    cols_del=[c for c in df.columns if 'Unnamed' in c]
    return df.drop(cols_del,axis=1)

def get_codon_table(aa, tax_id=None):
    """
    Gets host specific codon table.

    :param aa: list of amino acids
    :param host: name of host
    :returns: codon table (pandas dataframe)
    """
    
    # get codon table
    if tax_id is None:
        codontable=Data.CodonTable.unambiguous_dna_by_name["Standard"]
    else:
        codontable=Data.CodonTable.unambiguous_dna_by_id[tax_id]

    dcodontable=pd.DataFrame(pd.Series(codontable.forward_table))

    dcodontable.index.name='codon'
    dcodontable.columns=['amino acid']

    for cdn in codontable.stop_codons:
        dcodontable.loc[cdn,'amino acid']='*'        

    dcodontable=dcodontable.reset_index()
    rows=[]
    if isinstance(aa,list):
        for s in dcodontable['amino acid'].tolist():
            if s in aa:
                rows.append(True)
            else:
                rows.append(False)
    else:
        rows=dcodontable['amino acid']==aa
#     print(sum(rows))
    dcodontable=dcodontable.loc[rows,:].set_index('codon').reset_index()
    return dcodontable

def get_codon_usage(cuspp):
    """
    Creates codon usage table.

    :param cuspp: path to cusp generated file
    :returns: codon usage table (pandas dataframe)
    """
    # get codon usage stats
    dcodonusage=pd.read_csv(cuspp,sep='\t',header=5,keep_default_na=False)
    cols=''.join(dcodonusage.columns.tolist()).split(' ')
    dcodonusage.columns=[cols[-1]]
    dcodonusage.index.names=cols[:-1]

    dcodonusage=dcodonusage.reset_index().set_index('Codon')
    dcodonusage['amino acid']=[SeqUtils.seq1(s) for s in dcodonusage['#AA']]
    return dcodonusage

def saveemptytable(step,doutp=None):
    dout=pd.DataFrame(columns=stepi2cols[step])
    logging.warning(f"saved enmpty table {doutp}")
    if not doutp is None:
        dout.to_csv(doutp,sep='\t')
    else:
        return dout 

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

def get_possible_mutagenesis(dcodontable, dcodonusage, BEs, pos_muts):
    """
    Generates mutagenesis strategies from identities of reference and mutated codons (from dseq).
    :param cfg: configurations from yml file  
    """
    def write_dmutagenesis(cdni,posi,codon,codonmut,ntwt,ntmut,aa,aamut,method):
        """
        Write dmutagenesis table for each iteraction in get_possible_mutagenesis.
        """
        dmutagenesis.loc[cdni,'codon']=codon
        dmutagenesis.loc[cdni,'position of mutation in codon']=int(posi)
        dmutagenesis.loc[cdni,'codon mutation']=codonmut
        dmutagenesis.loc[cdni,'nucleotide']=ntwt
        dmutagenesis.loc[cdni,'nucleotide mutation']=ntmut
        dmutagenesis.loc[cdni,'amino acid']=aa
        dmutagenesis.loc[cdni,'amino acid mutation']=aamut
        dmutagenesis.loc[cdni,'mutation on strand']=method.split(' on ')[1]
        dmutagenesis.loc[cdni,'strand: mutation']=method.split(' on ')[1].replace(' strand','')
        dmutagenesis.loc[cdni,'method']=method.split(' on ')[0]                        
        dmutagenesis.loc[cdni,'codon mutation usage Fraction']=dcodonusage.loc[codonmut,'Fraction']
        dmutagenesis.loc[cdni,'codon mutation usage Frequency']=dcodonusage.loc[codonmut,'Frequency']
        return dmutagenesis
    
    def get_sm(dmutagenesis,BEs,positions,codon,muti,cdni):
        """
        Fetches single nucleotide mutagenesis strategies.
        """
        for method in BEs:
            for posi in positions: 
                if BEs[method][0]==codon[posi]:
                    for ntmut in BEs[method][1]:
                        if posi==0:
                            codonmut='{}{}{}'.format(ntmut,codon[1],codon[2])
                        elif posi==1:
                            codonmut='{}{}{}'.format(codon[0],ntmut,codon[2])
                        elif posi==2:
                            codonmut='{}{}{}'.format(codon[0],codon[1],ntmut)
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=1))
                        # if (aamut!='*') and (aamut!=aa): #  nonsence and synonymous
                        if muti==0:
                            cdni=cdni
                        else:
                            cdni=len(dmutagenesis)+1
                        muti+=1
                        ntwt=BEs[method][0]
                        if '-' in method.split(' on ')[1]:
                            ntwt=str(Seq.Seq(ntwt,Alphabet.generic_dna).reverse_complement())
                            ntmut=str(Seq.Seq(ntmut,Alphabet.generic_dna).reverse_complement())
                        dmutagenesis_row={'cdni':cdni,
                        'posi':posi+1,     #核苷酸在密码子中的位置
                        'codon':codon,      #原始密码子
                        'codonmut':codonmut,     #突变后的密码子
                        'ntwt':ntwt,     #野生型核苷酸
                        'ntmut':ntmut,    #突变后的核苷酸
                        'aa':aa,          #原始氨基酸
                        'aamut':aamut,     #突变后的氨基酸
                        'method':method}      #碱基编辑器编辑的方法
#                         print(dmutagenesis_row)
                        dmutagenesis=write_dmutagenesis(**dmutagenesis_row)
#                 else:
#                     logging.warning(f"BEs[{method}][0]!=codon[{posi}]")
        return dmutagenesis,muti
    
    #double nucleotide mutations
    positions={0:'@1st position',1:'@2nd position',2:'@3rd position'}
    #double nucleotide mutations
    positions_dm=[(i,j)  for i in positions.keys() for j in positions.keys() if i<j]
    #double nucleotide mutations
    positions_tm=[[0,1,2]]

    dmutagenesis=dcodontable.copy()

    for cdni in dmutagenesis.index:
            codon=dmutagenesis.loc[cdni,'codon']
            aa=dmutagenesis.loc[cdni,'amino acid']
            muti=0
            #single nucleuotide mutations
            dmutagenesis,muti=get_sm(dmutagenesis,BEs,positions,codon,muti,cdni)
    if len(dmutagenesis)==0:
        logging.warning('no guides designed; saving an empty table.')
        dmutagenesis=saveemptytable(step=2)
    else:
        dmutagenesis=dmutagenesis.dropna()  
        dmutagenesis['nucleotide mutation: count']=[len(s) for s in dmutagenesis['nucleotide mutation']]
        dmutagenesis=dmutagenesis.sort_values('codon') 
        dmutagenesis=dmutagenesis.set_index('method').join(pos_muts)
        dmutagenesis=dmutagenesis.reset_index()

        nt2complement=get_nt2complement()
        dmutagenesis['nucleotide: wild-type']=dmutagenesis.apply(lambda x : x['nucleotide'] if x['strand: mutation']=='+' else reverse_complement_multintseq(x['nucleotide'],nt2complement),axis=1) 
        dmutagenesis['nucleotide: mutation']=dmutagenesis.apply(lambda x : x['nucleotide mutation'] if x['strand: mutation']=='+' else reverse_complement_multintseq(x['nucleotide mutation'],nt2complement),axis=1)
    
    return dmutagenesis

def filterdmutagenesis(dmutagenesis, mutation_type, keep_mutation_nonsense, max_subs_per_codon, BE_names, mutations, non_intermutables):
    """
    Filters the mutagenesis strategies by multiple options provided in configuration file (.yml).

    :param dmutagenesis: mutagenesis strategies (pd.DataFrame)
    :param cfg: configurations from yml file
    """
    logging.info('filtering: dmutagenesis.shape: '+str(dmutagenesis.shape))
    
    # filter by mutation_type
    if not mutation_type is None:   #是否是允许同义突变
        if mutation_type=='S':       
            dmutagenesis=dmutagenesis.loc[(dmutagenesis['amino acid']==dmutagenesis['amino acid mutation'])]
        elif mutation_type=='N':
            dmutagenesis=dmutagenesis.loc[(dmutagenesis['amino acid']!=dmutagenesis['amino acid mutation'])]
    logging.info('filtering by mutation_type: dmutagenesis.shape: '+str(dmutagenesis.shape))
    
    # filter by nonsense
    if not keep_mutation_nonsense is None:
        if not keep_mutation_nonsense:
            dmutagenesis=dmutagenesis.loc[(dmutagenesis['amino acid mutation']!='*'),:]
    logging.info('filtering by nonsense: dmutagenesis.shape: '+str(dmutagenesis.shape))
    
    # filter by mutation per codon
    if not max_subs_per_codon is None:
        dmutagenesis=dmutagenesis.loc[(dmutagenesis['nucleotide mutation: count']==max_subs_per_codon ),:]
    logging.info('filtering by mutation per codon: dmutagenesis.shape: '+str(dmutagenesis.shape))
    
    # filter by method
    if not BE_names is None:
        dmutagenesis=dmutagenesis.loc[dmutagenesis['method'].isin([BE_names]),:]
        logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    
    logging.info('filtering by method: dmutagenesis.shape: '+str(dmutagenesis.shape))
    
    # filter by submap
    if (mutations=='mimetic') or (mutations=='substitutions'):
            if mutations == 'mimetic':                
                dsubmap=get_submap_mimetic(cfg)
            elif mutations == 'substitutions':
                dsubmap=pd.read_csv(cfg['dsubmap_preferred_path'],sep='\t',keep_default_na=False) # has two cols: amino acid and amino acid mutation
            import seaborn as sns
            dsubmap.to_csv(f"{cfg['datad']}/dsubmap.tsv",sep='\t')
            dmutagenesis=pd.merge(dsubmap,dmutagenesis,on=['amino acid','amino acid mutation'],how='inner')
            dsubmap['count']=1
            sns.heatmap(dsubmap.pivot_table(columns='amino acid',index='amino acid mutation',values='count'),square=True)
            plt.xlabel('wild-type amino acid')
            plt.ylabel('mutated amino acid')
            plt.savefig(f"{cfg['datad']}/heatmap_submap.svg")
            logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    
    logging.info('filtering by submap: dmutagenesis.shape: '+str(dmutagenesis.shape))
    
    # filter by non interchageables
    if not non_intermutables is None:
        if len(non_intermutables)!=0:               
            dmutagenesis=dmutagenesis.loc[~(dmutagenesis['amino acid'].isin(cfg['non_intermutables']) \
                                               & dmutagenesis['amino acid mutation'].isin(cfg['non_intermutables'])),:]
            logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    
    logging.info('filtering by non interchageables: dmutagenesis.shape: '+str(dmutagenesis.shape))
    return dmutagenesis



def get_mutations(dsequencesp, dbepamsp, dmutagenesisp, dmutagenesisallp):
   
    if not exists(dmutagenesisp):
        dseq=del_Unnamed(pd.read_csv(dsequencesp,sep='\t',keep_default_na=False))
        dcodontable=get_codon_table(aa=aminoacids)
        dcodonusage=get_codon_usage(cuspp=cuspp)
        
        dbepams=del_Unnamed(pd.read_table(dbepamsp,keep_default_na=False))
        
        dBEs=dbepams.loc[:,cols_dbes]
        print(dBEs['method'].unique())
        BEs2mutations={}
        for method in dBEs['method'].unique():
            for strand in dBEs['strand'].unique():
                dBEsi=dBEs.loc[(dBEs['method']==method) & (dBEs['strand']==strand),:]
                BEs2mutations[f"{method} on {strand} strand"]=[dBEsi['nucleotide'].unique().tolist()[0],
                                                        dBEsi['nucleotide mutation'].unique().tolist()]
                
        pos_muts=dBEs.loc[:,['method']+['distance of mutation from PAM: minimum',
                'distance of mutation from PAM: maximum',
                'distance of codon start from PAM: minimum',
                'distance of codon start from PAM: maximum']].drop_duplicates().set_index('method')

        dmutagenesis = get_possible_mutagenesis(dcodontable=dcodontable, dcodonusage=dcodonusage, BEs=BEs2mutations, pos_muts=pos_muts)
        dmutagenesis.to_csv(dmutagenesisallp,sep='\t')


        mutation_type = 'both'
        keep_mutation_nonsense = False
        max_subs_per_codon = 1
        BE_names = 'Target-AID'
        mutations = 'mutations'


        dmutagenesis=filterdmutagenesis(dmutagenesis, mutation_type, keep_mutation_nonsense, max_subs_per_codon, BE_names, mutations, non_intermutables=None)
        colns_pos=[c for c in dmutagenesis if ('position' in c) or ('Position' in c)]
                
                
        dmutagenesis.loc[:,colns_pos].astype(int)
        dmutagenesis.to_csv(dmutagenesisp,sep='\t')


        cols_dpam=['PAM', 'PAM position', 'guide length']
        dmutagenesis=dmutagenesis.merge(dbepams.loc[:,['method']+cols_dpam],on='method')

        logging.info('Possible 1 nucleotide mutations:')
        logging.info(dmutagenesis.set_index('amino acid')[['amino acid mutation','method','codon','codon mutation',
                #               'position of mutation in codon','mutation on strand',
                #               'nucleotide','nucleotide mutation',
                            ]])
        for aa in aminoacids:
            logging.info(aa+' can be mutated to:')
            logging.info(list(dmutagenesis.loc[dmutagenesis.loc[:,'amino acid']==aa,:].loc[:,'amino acid mutation'].unique()))
                
    import gc
    gc.collect()  



def main():
    

    p_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))

    basedir = os.path.join(p_dir, 'data/YB01/02_mutagenesis/')

    if not exists(basedir):
        os.makedirs(basedir)


    dsequencesp = os.path.join(p_dir,'data/YB01/01_sequences/dsequences.tsv')
    dbepamsp = '/home/yanghe/githubcode/beditor/sgRNA_primer/input_output/YB01/00_input//dbepams.tsv'

    dmutagenesisp = os.path.join(basedir,'dmutagenesis.tsv')
    dmutagenesisallp = os.path.join(basedir,'dmutagenesis_all.tsv')

    get_mutations(dsequencesp, dbepamsp, dmutagenesisp, dmutagenesisallp) 


   


if __name__ == '__main__':

    main()