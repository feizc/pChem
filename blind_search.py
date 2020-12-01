# 在人工筛选常见修饰后，对数据进行盲搜鉴定

import os 
from utils import parameter_file_read, modification_ini_path, modification_ini_dict, \
    modification_ini_generation, blind_cfg_write, search_exe_path, mass_diff_list_generate


current_path = os.getcwd()
pchem_cfg_path = os.path.join(os.getcwd(), 'pChem.cfg')
blind_cfg_path = os.path.join(os.path.join(os.getcwd(), 'template'), 'blind.cfg')

if __name__ == "__main__":
    parameter_dict = parameter_file_read(pchem_cfg_path)

    # 读取所有的modification
    modification_path = modification_ini_path(parameter_dict)
    modification_dict = modification_ini_dict(modification_path)

    # 重新生成modification.ini文件
    common_modification_list = modification_ini_generation(current_path, modification_dict)

    # 重新生成blind.cfg文件
    res_path = blind_cfg_write(blind_cfg_path, current_path, parameter_dict, common_modification_list)
    '''
    # 运行盲搜search.exe进行搜索
    bin_path, exe_path = search_exe_path(parameter_dict)
    cmd = exe_path + ' ' + blind_cfg_path 
    os.chdir(bin_path)
    receive = os.system(cmd)
    print(receive)
    '''

    # 读取鉴定结果，生成位置修饰的候选列表
    mass_diff_list_generate(res_path, current_path)
