# 限定式搜索
# 1. 对未知修饰的质量做校正
# 2. 运行限定式搜索输出最终鉴定结果

import os 
from utils import parameter_file_read, mass_diff_read, expand_modification_ini, modification_ini_path, close_cfg_write, search_exe_path
from mass_diff_correction import mass_correct, small_delta_filter, mass_diff_diff_filter  


def close_search():
    
    current_path = os.getcwd()
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    close_cfg_path = os.path.join(os.path.join(current_path, 'template'), 'close.cfg')

    # 读取位置修饰质量数的列表
    mass_diff_list = mass_diff_read(current_path)

    
    # 盲搜的结果文件
    parameter_dict = parameter_file_read(pchem_cfg_path)
    blind_path = os.path.join(parameter_dict['output_path'], 'blind')
    blind_path = os.path.join(blind_path, 'pFind-Filtered.spectra')
    
    name2mass, mass_diff_list = small_delta_filter(mass_diff_list)
    refined_list = mass_diff_diff_filter(name2mass, mass_diff_list, parameter_dict['mass_diff_diff'])
    print(refined_list)

    '''
    # 对未知质量数质量做校正
    mass_diff_dict = mass_correct(current_path, blind_path, mass_diff_list)

    # 统计未知质量发生的位点
    # 这里暂时假设发生在N端 'PEP_N' or 'NORMAL'
    position_dict = {}
    for key in mass_diff_dict.keys():
        position_dict[key] = 'PEP_N'

    # 修改的modification-new.ini文件，加入质量数修饰
    ini_path = modification_ini_path(parameter_dict)
    new_ini_path = expand_modification_ini(mass_diff_dict, position_dict, current_path, ini_path)

    # 生成限定式参数文件 
    res_path = close_cfg_write(close_cfg_path, current_path, parameter_dict, mass_diff_dict)
    # 调用pfind进行限定式搜索 
    bin_path, exe_path = search_exe_path(parameter_dict)
    cmd = exe_path + ' ' + close_cfg_path 
    os.chdir(bin_path)
    receive = os.system(cmd)
    print(receive)
    '''

if __name__ == "__main__":
    close_search()
