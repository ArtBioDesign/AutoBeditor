import sys,os
# sys.path.append('../')  
from primerDesign.src.sgRNA_utils import sgRNA_primer_util  as util
from primerDesign.src.sgRNA_utils.sgRNA_primer_config import config



#生成测序引物模板
def design_sequencing_primers(dict_plasmid_id, dict_plasmid_seq,mute_position=0):
    '''
        生成测序引物模板（在target_gene_seq 前后加200bp）,设计测序引物
        params:     
            dict_plasmid_id  
            dict_plasmid_seq
        returns:
    '''
    #生成引物模板
    target_gene_down_seq = dict_plasmid_seq['target_gene_down_seq']
    target_gene_up_seq = dict_plasmid_seq['target_gene_up_seq']
    target_gene_seq = dict_plasmid_seq['mute_after_target_gene_seq']

    
    # target_gene_seq 前后加200bp    
    if len(target_gene_down_seq) >= 200:
        sequencing_peimers_template = target_gene_seq + target_gene_down_seq[:200]
    else:
        sequencing_peimers_template = target_gene_seq + target_gene_down_seq
        temp_len = 200 - len(target_gene_down_seq) 
        sequencing_peimers_template = sequencing_peimers_template + target_gene_up_seq[:temp_len]

    if len(target_gene_up_seq) >= 200:
        sequencing_peimers_template = target_gene_up_seq[-200:] + sequencing_peimers_template
    else:
        sequencing_peimers_template = target_gene_up_seq + sequencing_peimers_template
        temp_len = 200 - len(target_gene_up_seq)
        sequencing_peimers_template = target_gene_down_seq[-temp_len:] + sequencing_peimers_template
    #设计引物 
    if mute_position ==0:
        result = design_seq_primer(seq_id=dict_plasmid_id, seq_temp=sequencing_peimers_template, seq_target=target_gene_seq)
        print(result)
        
    else:

        result = design_seq_primer(seq_id=dict_plasmid_id, seq_temp=sequencing_peimers_template, seq_target=target_gene_seq,mute_position=mute_position)
         
    return result

#设计测序引物
def design_seq_primer(seq_id, seq_temp, seq_target, mute_position=0):
    dict_seq_primer={}
    dict_seq_primer_failtrue={}
    len_target=len(seq_target)     
    site_target_temp=seq_temp.find(seq_target)
    temp_p1=seq_temp[site_target_temp-120:site_target_temp-80]

    dict_res_p1 = util.primer_design(seqId=seq_id, 
                                        seqTemplate=temp_p1,
                                        stype='none',
                                        mute_type='sequencing',
                                        global_args=config.Q_GLOBAL_ARGS
                                    )
   
    #第一条引物
    judge_primer_is_or_not( seq_id,dict_res_p1,
                            dict_seq_primer,
                            dict_seq_primer_failtrue,
                            primer_name="SEQUENCING_PRIMER_1",
                            type='LEFT')

    if 600 < len_target <= 1200:
        temp_p2 = seq_temp[site_target_temp + len_target + 80 : site_target_temp + len_target + 120]
        dict_res_p2=util.primer_design(seqId=seq_id,
                                        seqTemplate=temp_p2,
                                        stype='none',
                                        mute_type='sequencing',
                                        global_args=config.Q_GLOBAL_ARGS
                                    )
         #判断引物是否成功
        judge_primer_is_or_not( seq_id,dict_res_p2,
                            dict_seq_primer,
                            dict_seq_primer_failtrue,
                            primer_name="SEQUENCING_PRIMER_2",
                            type='RIGHT')

    elif 1200 < len_target:    
 
        temp_p2 = seq_temp[site_target_temp+len_target+80:site_target_temp+len_target+120]
        dict_res_p2=util.primer_design(seqId=seq_id,
                                        seqTemplate=temp_p2,
                                        stype='none',
                                        mute_type='sequencing',
                                        global_args=config.Q_GLOBAL_ARGS
                                    )  
         #判断引物是否成功
        judge_primer_is_or_not( seq_id,dict_res_p2,
                                dict_seq_primer,
                                dict_seq_primer_failtrue,
                                primer_name="SEQUENCING_PRIMER_2",
                                type='RIGHT')
    
        k=2
        while True:
            distance = [int(site_target_temp)+480, int(site_target_temp)+520]
            if mute_position !=0:
                # distance = [int(site_target_temp)+480, int(site_target_temp)+520]
                noisy_data = [i for i in mute_position if  (i+200) >=distance[0] and (i+200) <=distance[1]]

                init_position1 =  distance[0]    
                init_position2 =  distance[1]
                times = 0
                while len(noisy_data) != 0:
                    distance[0] = distance[0] - 1
                    distance[1] = distance[1] - 1
                    noisy_data = [i for i in mute_position if  (i+200) >=distance[0] and (i+200) <=distance[1]]
                    times = times + 1 
                
                    if times == 199:
                        distance[0] = init_position1
                        distance[1] = init_position2  
                        break
    
                temp_pn = seq_temp[distance[0]: distance[1]]
            else:
                temp_pn = seq_temp[distance[0]: distance[1]]
            dict_res_pn = util.primer_design(seqId=seq_id,
                                                seqTemplate=temp_pn,
                                                stype='none',
                                                mute_type='sequencing',
                                                global_args=config.Q_GLOBAL_ARGS
                                            )


            #上一次测序长度 = 0 ----distance[1]+80
            last_sequencing_seq = seq_temp[site_target_temp:distance[1]+80]
            last_sequencing_len = distance[1]+80 - site_target_temp
            #确定此次引物的作用域的起始到最后
            seq_target = seq_target[distance[1]-200+80:]  
            site_target_temp=seq_temp.find(seq_target) 
            len_target=len(seq_target)

            if len_target>=600:
                k+=1
                     #判断引物是否成功
                judge_primer_is_or_not( seq_id,dict_res_pn,
                                        dict_seq_primer,
                                        dict_seq_primer_failtrue,
                                        primer_name=f"SEQUENCING_PRIMER_{k}",
                                        type='LEFT')
            else:
                break

        
        if 200 <= len_target < 600:
            length_add=int((600-len_target)/2)
            temp_pe=seq_temp[site_target_temp-length_add-120:site_target_temp-length_add-80]
            dict_res_pe = util.primer_design(seqId=seq_id,
                                            seqTemplate=temp_pe,
                                            stype='none',
                                            mute_type='sequencing',
                                            global_args=config.Q_GLOBAL_ARGS
                                        )
            #判断引物是否成功
            judge_primer_is_or_not( seq_id,dict_res_pe,
                                    dict_seq_primer,   
                                    dict_seq_primer_failtrue,
                                    primer_name=f"SEQUENCING_PRIMER_{k}",     
                                    type='LEFT')

    return dict_seq_primer, dict_seq_primer_failtrue

#判断引物是与否
def judge_primer_is_or_not( seq_id, dict_res,primer_suc,primer_fail,primer_name,type='LEFT'):
    if len(list(dict_res.keys())) < 10:
        primer_fail[f'{seq_id}:{primer_name}'] = dict_res   
        # print(seq_id,primer_name, dict_res) 
        primer_suc[primer_name] = {'failtrue':dict_res}   
    else:
        primer_suc[primer_name] = dict_res[f'PRIMER_{type}_0_SEQUENCE']
        primer_suc[f'{primer_name}_TM'] =dict_res[f'PRIMER_{type}_0_TM']


              
  
    

