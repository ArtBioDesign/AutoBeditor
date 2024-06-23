
from glob import glob
import os
import numpy as np
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
import logging
from os import makedirs
import pandas as pd 

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

def get_genomes_by_path(host,genome_release,genome_fna_path,genome_gff_path):
    '''
    cd -f genome_fna_path,genome_gff_path  genomep,genomegffp
    bwa index genomep
    samtools faidx genomep
    cut -f1,2 genomep.fai > genomep.sizes

    '''

    host_="_".join(s for s in host.split('_')).capitalize()

    ncbi_fastad = 'pub/{}/fasta/{}/dna/'.format(genome_release, host)
    genome_fastad = '{}/{}'.format(dirname(realpath(__file__)),ncbi_fastad)
    genomep = '{}genome.fna'.format(genome_fastad)
    if not exists(genomep): 
        if not exists(genome_fastad):
            os.makedirs(f'{genome_fastad}')
        gene =  '{}{}'.format(genome_fastad,'genome.fna')
        cmd=f"cp -f {genome_fna_path} {gene}"
        os.system(cmd)


    if not exists(genomep+'.bwt'):
        cmd='{} index {}'.format('bwa',genomep)
        print(cmd)
        os.system(cmd)
    else:        
        logging.info('bwa index is present')
    if not exists(genomep +'.fai'):
        cmd='{} faidx {}'.format('samtools',genomep)
        print(cmd)
        os.system(cmd)
    else:
        logging.info('samtools index is present')
    if not exists(genomep+'.sizes'):
        cmd='cut -f1,2 {}.fai > {}.sizes'.format(genomep,genomep)
        print(cmd)
        os.system(cmd)
    else:
        logging.info('sizes of contigs are present')

    #gff
    ncbi_gff3d = 'pub/{}/gff3/{}/'.format(genome_release,host)                            
    genome_gff3d=f'{dirname(realpath(__file__))}/{ncbi_gff3d}'
    genomegffp='{}/genome.gff'.format(genome_gff3d)
    if not exists(genomegffp):
        if not exists(genome_gff3d):
            os.makedirs(f'{genome_gff3d}')
        cmd=f"cp -f {genome_gff_path} {genomegffp}"
        os.system(cmd)

    return genomep, genomegffp

def dbes2dbes_strands(dBEs):
    """
    add reverse strand editing methods in the dBEs dataframe   添加反义链的编辑方法
    
    :param dBEs: pandas dataframe with BE methods       
    """

    dBEs_=dBEs.copy()
    dBEs_['strand']='-'
    nt2complement=get_nt2complement()
    for col_nt in ['nucleotide','nucleotide mutation']:
        dBEs_[col_nt]=dBEs_[col_nt].apply(lambda x : nt2complement[x])
    dBEs=dBEs.append(dBEs_,sort=True)
    return dBEs

def to_table(df,p):
    if p.endswith('.tsv') or p.endswith('.tab'):
        if not exists(dirname(p)) and dirname(p)!='':
            makedirs(dirname(p),exist_ok=True)
        df.to_csv(p,sep='\t')
    elif p.endswith('.pqt') or p.endswith('.parquet'):
        to_table_pqt(df,p)
    else: 
        logging.error(f'unknown extension {p}')   
        
def to_table_pqt(df,p):
    if not exists(dirname(p)) and dirname(p)!='':
        makedirs(dirname(p),exist_ok=True)
    df.to_parquet(p,engine='fastparquet',compression='gzip',)

def run( dinp, genome_fna_path, genome_gff_path, cfpg, param):

    #0    
    genome_release = param.get("genome_release")
    host = param.get("host")
    be_names = param.get("be_names")
    pams = param.get("pams")


    #1
    genomep, genomegffp = get_genomes_by_path(host,genome_release,genome_fna_path,genome_gff_path)

    #2
    prjd=os.path.join(cfpg,host)
    steps = ['00_input', '01_sequence', '02_mutagenesis', '03_guides', '04_offtargets']
    for i in range(5):
        step = os.path.join(prjd,steps[i])
        steps[i] = step
        if not exists(step):
            makedirs(step)
    
    #3
    dbepamsp = os.path.join(steps[0],'dbepams.tsv')
    dbepams = pd.read_table(f"{dirname(realpath(__file__))}/input/data/dbepams.tsv",keep_default_na=False)
    dbepams['strand']='+'
    dbepams=dbes2dbes_strands(dbepams)
    dbepams=dbepams.loc[(dbepams['method'].isin(be_names) & dbepams['PAM'].isin(pams) ),:]
    to_table(dbepams,dbepamsp)

    #4
    flankntc = 22
    bedp = os.path.join(steps[1], 'dbedntmuts.bed')
    fastap = os.path.join(steps[1], 'dbedntmuts.fa')
    dbedntmutsp = os.path.join(steps[1], 'dbedntmuts.tsv')
    dsequencesp = os.path.join(steps[1], 'dsequences.tsv')
    from beditor.lib.get_seq_nucletide import get_seq_nucleotide
    dsequences = get_seq_nucleotide(genomep, dinp, bedp, fastap, dbedntmutsp, dsequencesp,flankntc)

    #5
    dsequencesp = os.path.join(steps[1],'dsequences.tsv')  
    dmutagenesisp = os.path.join(steps[2],'dmutagenesis.tsv')
    dmutagenesisallp = os.path.join(steps[2],'dmutagenesis_all.tsv')
    from beditor.lib.get_mutations import get_mutations
    get_mutations(dsequencesp, dbepamsp, dmutagenesisp, dmutagenesisallp)

    #6
    dguideslinp = os.path.join(steps[3], 'dguides.tsv')      
    dguides_nofltp = os.path.join(steps[3], 'dguides_noflt.tsv')
    dmutagenesisp_03 = os.path.join(steps[3], 'dmutagenesis.tsv')  
    dmutagenesisp_02 = os.path.join(steps[2], 'dmutagenesis.tsv')
    dsequences_01 = os.path.join(steps[1], 'dsequences.tsv')
    dsequences_03 = os.path.join(steps[3], 'dsequences.tsv')
    mutation_format = 'nucleotide'    
    flankaas=7
    make_control_pos = None
    make_control_neg = None
    from beditor.lib.get_guides import get_guides
    get_guides( dguideslinp, dguides_nofltp, dmutagenesisp_03, dmutagenesisp_02, dsequences_01, dsequences_03, mutation_format, dbepamsp, flankaas, make_control_pos, make_control_neg )

    #7
    dofftargetsp = os.path.join(steps[4], 'dofftargets.tsv') 
    dguidesp = os.path.join(steps[3],'dguides.tsv')
    datatmpd = os.path.join(steps[4], 'tmp/')
    if not exists(datatmpd):
        os.makedirs(datatmpd)
    from beditor.lib.get_dguides2offtargets import get_dguides2offtargets
    get_dguides2offtargets(genomep, genomegffp, datatmpd, dguidesp, dofftargetsp, dbepamsp)



def main():

    dinp = './input/YB01/YB01_beditor_input.tsv'
    genome_fna_path = './input/YB01/YB01.fna'  
    genome_gff_path = './input/YB01/YB01_v2.gff'
    workdir = './output'

    param={ 
        "genome_release" : 'GCF001',
        "host" : 'yb01',
        "be_names" : ['Target-AID'],
        "pams" : ['NG']
    }
    
    run(dinp, genome_fna_path, genome_gff_path, workdir, param )


if __name__ == '__main__':

    main()


