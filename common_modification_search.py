# 使用开放式寻找常见的修饰列表
import os
from utils import parameter_file_read, open_cfg_write, search_exe_path, spectra_result_read

current_path = os.getcwd()
pchem_cfg_path = os.path.join(os.getcwd(), 'pChem.cfg')
open_cfg_path = os.path.join(os.path.join(os.getcwd(), 'template'), 'open.cfg')



if __name__ == "__main__":
    # 读取参数
    parameter_dict = parameter_file_read(pchem_cfg_path)
    
    # 写入开放式pFind参数文件
    res_path = open_cfg_write(open_cfg_path, parameter_dict, 'open')
'''
    # 调用pFind进行搜索
    bin_path, exe_path = search_exe_path(parameter_dict)
    cmd = exe_path + ' ' + open_cfg_path 
    os.chdir(bin_path)
    receive = os.system(cmd)
    print(receive)
    
    # 读取结果文件，给出常见修饰列表
    spectra_res_path = os.path.join(res_path, 'pFind.summary')
    spectra_result_read(spectra_res_path, current_path)
'''
