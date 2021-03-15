# 中性丢失鉴定 

from utils import parameter_file_read 
import os 

def ion_type_determine(current_path): 
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    parameter_dict = parameter_file_read(pchem_cfg_path) 
    print(parameter_dict)


if __name__ == "__main__": 
    current_path = os.getcwd() 
    ion_type_determine(current_path)