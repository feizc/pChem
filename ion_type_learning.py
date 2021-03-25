# 中性丢失鉴定 

from utils import parameter_file_read 
import os 
from mass_diff_correction import amino_acid_dict, common_dict_create 

# 存放谱图的类
class MassSpectrum:
    def __init__(self, charge, pepmass, peak_list):
        self.charge = charge 
        self.pepmass = pepmass
        self.peak_list = peak_list 


# 从mgf文件中读取谱图
def mgf_read(mgf_path):
    mass_spectra_dict = {}
    with open(mgf_path, 'r') as f:
        lines = f.readlines()
    # print(len(lines))
    i = 0 
    while i < len(lines):
        if 'BEGIN' in lines[i]: 
            i += 1
            spectrum_name = lines[i].split('=')[1].strip()
            i += 1
            spectrum_charge = int(lines[i].split('=')[1][0])
            i += 2 
            spectrum_pepmass = float(lines[i].split('=')[1])
            spectrum_peak_list = []
            while i < len(lines):
                i += 1 
                if 'END' in lines[i]:
                    break 
                spectrum_peak_list.append(lines[i].split()) 
            # print(spectrum_peak_list)
        spectrum = MassSpectrum(spectrum_charge, spectrum_pepmass, spectrum_peak_list)
        mass_spectra_dict[spectrum_name] = spectrum 
        i += 1 
        break 
    # print('The number of spectra: ', len(mass_spectra_dict.keys()))
    return mass_spectra_dict


# 读取盲搜的结果 
def blind_res_read(blind_res_path):
    # print(blind_res_path) 
    with open(blind_res_path, 'r') as f:
        lines = f.readlines() 
    return lines[1:]


# 筛选出含有指定修饰的PSM
def PSM_filter(blind_res, modification): 
    filtered_res = []
    for line in blind_res: 
        if modification in line:
            filtered_res.append(line) 
    return filtered_res 


# 统计中性丢失的数目 
def ion_type_compute(filtered_res, modification):
    print(filtered_res[0])



# 离子类型学习 
def ion_type_determine(current_path, modification_list): 
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    parameter_dict = parameter_file_read(pchem_cfg_path) 
    # print(parameter_dict)

    # 质谱数据读取 
    mass_spectra_dict = {} 
    for msms_path in parameter_dict['msms_path']:
        mgf_path = msms_path.split('=')[1].split('.')[0] + '.mgf'
        cur_mass_spectra_dict = mgf_read(mgf_path) 
        mass_spectra_dict.update(cur_mass_spectra_dict)
    print('The number of spectra: ', len(mass_spectra_dict.keys())) 

    # 读取盲搜得到的结果 
    blind_path = os.path.join(parameter_dict['output_path'], 'blind')
    blind_res_path = os.path.join(blind_path, 'pFind-Filtered.spectra')
    blind_res = blind_res_read(blind_res_path) 

    # 读取常见修饰的列表 
    common_modification_dict = common_dict_create(current_path)
    # print(common_modification_dict)

    # 筛选有效的PSM 
    for modification in modification_list:
        filtered_res = PSM_filter(blind_res, modification) 
        ion_type_compute(filtered_res, modification)
        # print(filtered_res)






if __name__ == "__main__": 
    current_path = os.getcwd() 
    modification_list = ['PFIND_DELTA_252', 'PFIND_DELTA_258']
    modification_dict = {'PFIND_DELTA_252': 252.121858, 'PFIND_DELTA_258':258.141955}
    ion_type_determine(current_path, modification_list)
    