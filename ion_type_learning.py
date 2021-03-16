# 中性丢失鉴定 

from utils import parameter_file_read 
import os 

class MassSpectrum:
    def __init__(self, charge, pepmass, peak_list):
        self.charge = charge 
        self.pepmass = pepmass
        self.peak_list = peak_list 


def ion_type_determine(current_path): 
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    parameter_dict = parameter_file_read(pchem_cfg_path) 
    
    mass_spectra_dict = {}
    # 数据读取 
    for msms_path in parameter_dict['msms_path']:
        print(msms_path)
    


# 从mgf文件中读取谱图
def mgf_read():
    print('1')


if __name__ == "__main__": 
    current_path = os.getcwd() 
    ion_type_determine(current_path)