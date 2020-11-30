import os 

# 读取参数行内容并返回
def parameter_pick(line):
    eq_idx = line.find('=')
    parameter_content = line[eq_idx+1:].strip()
    return parameter_content

# 改写参数行的内容并返回
def parameter_modify(line, content):
    eq_idx = line.find('=')
    line = line[:eq_idx+1]
    line += content
    line += '\n'
    return line

# 读取参数文件
def parameter_file_read(path):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    parameter_dict = {}
    for line in lines:
        if 'pfind_install' in line:
            parameter_dict['pfind_install_path'] = parameter_pick(line)
        if 'fasta_path' in line:
            parameter_dict['fasta_path'] = parameter_pick(line)
        if 'msms_path' in line:
            parameter_dict['msms_path'] = parameter_pick(line)
        if 'output_path' in line:
            parameter_dict['output_path'] = parameter_pick(line)
            if not os.path.exists(parameter_dict['output_path']):
                os.mkdir(parameter_dict['output_path'])
    return parameter_dict

# 写参数文件
def pfind_cfg_write(cfg_path, parameter_dict, pattern):
    with open(cfg_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 指定参数内容修改
    for i in range(len(lines)):
        if 'modpath' in lines[i]:
            mod_path = os.path.join(parameter_dict['pfind_install_path'], 'bin')
            mod_path = os.path.join(mod_path, 'modification.ini')
            lines[i] = parameter_modify(lines[i], mod_path)
        if 'fastapath' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['fasta_path'])
        if 'outputpath' in lines[i]:
            res_path = os.path.join(parameter_dict['output_path'], pattern)
            if not os.path.exists(res_path):
                os.mkdir(res_path)
            lines[i] = parameter_modify(lines[i], res_path)
        if 'outputname' in lines[i]:
            lines[i] = parameter_modify(lines[i], 'open')
        if 'msmspath' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['msms_path'])
    
    # 写入参数文件
    with open(cfg_path, 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line)
    return res_path


# 返回search.exe地址
def search_exe_path(parameter_dict):
    bin_path = os.path.join(parameter_dict['pfind_install_path'], 'bin')
    exe_path = os.path.join(bin_path, 'Searcher.exe')
    return bin_path, exe_path 


# 读取谱图结果文件, 将常见的修饰列表写出来提供删选
def spectra_result_read(spectra_res_path, target_path):
    # 保存常见的修饰
    model_res = []
    with open(spectra_res_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if 'Modifications:' in lines[i]:
            i += 1
            model_res.append(lines[i])
            i += 1
            while True:
                freq = lines[i].split('\t')[1]
                idx_left = freq.find('(')
                idx_right = freq.find('%')
                freq = float(freq[idx_left+1:idx_right])
                if freq < 1.5:
                    break
                i += 1
                model_res.append(lines[i])
    # print(model_res)
    target_path = os.path.join(target_path, 'common_modification_list.txt')
    with open(target_path, 'w', encoding='utf-8') as f:
        for line in model_res:
            f.write(line)


