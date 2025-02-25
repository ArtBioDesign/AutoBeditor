import os

class config: 

        #单双点设计引物的全局参数  
    GLOBAL_ARGS = {
                'PRIMER_PICK_ANYWAY':1,
                'PRIMER_PRODUCT_SIZE_RANGE': 0,
                'PRIMER_NUM_RETURN':2
        }       

    #测序引物设计的全局参数
    Q_ARGS = {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_SIZE': 18,   
                'PRIMER_MAX_SIZE': 25,
                'PRIMER_OPT_TM': 65.0,
                'PRIMER_MIN_TM': 55.0,
                'PRIMER_MAX_TM': 75.0,    
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
    }

    Q_GLOBAL_ARGS = {   
                'PRIMER_PICK_ANYWAY':1,    
                'PRIMER_TASK': 'pick_primer_list', 
        }  
      



    INPUT_FILE_PATH = ''   
    OUTPUT_FILE_PATH =''    

    PLASMID_FILE_NAME = ''    
    INPUT_FILE_NAME = ''  
    NEW_OUTPUT_FILE_PATH=''
    
    
    SEQUENCING_PRIMER_SUCCESS = 'sequencing_primer_success.xlsx'
    SEQUENCING_PRIMER_FAILTRUE = 'sequencing_primer_failtrue.xlsx'

    OUTPUT_FILE_NAME_PLASMID_MUTATION = "plasmid_mutation.gb" 
    DATA_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) + '/data'
