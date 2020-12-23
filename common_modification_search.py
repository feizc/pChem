# 使用开放式寻找常见的修饰列表
import os
from utils import parameter_file_read, open_cfg_write, search_exe_path, spectra_result_read


# 将整个开放式搜索过程包装成一个函数
def open_search():

    # 参数文件的路径
    current_path = os.getcwd()
    pchem_cfg_path = os.path.join(os.getcwd(), 'pChem.cfg')
    open_cfg_path = os.path.join(os.path.join(os.getcwd(), 'template'), 'open.cfg')
    
    # 读取参数
    parameter_dict = parameter_file_read(pchem_cfg_path)
    print('parameter_dict: ', parameter_dict) 

    if parameter_dict['open_flag'] == 'True':

        # 写入开放式pFind参数文件
        res_path = open_cfg_write(open_cfg_path, parameter_dict)
        
        # 调用pFind进行搜索
        bin_path, exe_path = search_exe_path(parameter_dict)
        cmd = exe_path + ' ' + open_cfg_path 
        os.chdir(bin_path)
        receive = os.system(cmd)
        print(receive)
        
        # 读取结果文件，给出常见修饰列表
        spectra_res_path = os.path.join(res_path, 'pFind.summary')
        spectra_result_read(spectra_res_path, current_path)
        

if __name__ == "__main__":
    open_search()