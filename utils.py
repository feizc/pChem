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

# 写open参数文件
def open_cfg_write(cfg_path, parameter_dict, pattern):
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


# 写blind参数文件
def blind_cfg_write(cfg_path, current_path, parameter_dict, common_modification_list):
    with open(cfg_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    mod_line = ""
    for mod in common_modification_list:
        mod_line += (mod + ';')
    mod_line = mod_line[:-1] + '\n'

    # 指定参数内容修改
    for i in range(len(lines)):
        if 'selectmod' in lines[i]:
            lines[i] = parameter_modify(lines[i], mod_line)
        if 'modpath' in lines[i]:
            mod_path = os.path.join(current_path, 'modification-null.ini')
            lines[i] = parameter_modify(lines[i], mod_path)
        if 'fastapath' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['fasta_path'])
        if 'outputpath' in lines[i]:
            res_path = os.path.join(parameter_dict['output_path'], 'blind')
            if not os.path.exists(res_path):
                os.mkdir(res_path)
            lines[i] = parameter_modify(lines[i], res_path)
        if 'outputname' in lines[i]:
            lines[i] = parameter_modify(lines[i], 'blind')
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


# 返回modification.ini地址
def modification_ini_path(parameter_dict):
    bin_path = os.path.join(parameter_dict['pfind_install_path'], 'bin')
    modification_ini_path = os.path.join(bin_path, 'modification.ini')
    return modification_ini_path


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


# 读取所有候选修饰，并返回dict
def modification_ini_dict(path):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    modification_dict = {}
    i = 1
    while i < len(lines):
        # 防止文件后面的空行
        if len(lines) < 4:
            break
        modification_name = lines[i].split()[0]
        eq_idx = modification_name.find('=')
        modification_name = modification_name[eq_idx+1:]
        modification_dict[modification_name] = lines[i] + lines[i+1]
        i += 2
    return modification_dict


# 生成新的modification-null.ini用于检索
def modification_ini_generation(path, modification_dict):
    # 读取之前确定的常见修饰列表
    common_list_path = os.path.join(path, 'common_modification_list.txt')
    with open(common_list_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    lines = lines[1:]
    common_modification_list = []
    modification_ini_lines = []
    for line in lines:
        if len(line) < 4:
            break
        mod_name = line.split('\t')[0]
        common_modification_list.append(mod_name)
        modification_ini_lines.append(modification_dict[mod_name])
    # print(common_modification_list)
    # print(modification_ini_lines)
    
    # 写入新的modification-null.ini
    new_ini_path = os.path.join(path, 'modification-null.ini')
    i = 1
    with open(new_ini_path, 'w', encoding='utf-8') as f:
        f.write('@NUMBER_MODIFICATION=' + str(len(common_modification_list)) + '\n')
        for line in modification_ini_lines:
            if 'name' in line:
                eq_idx = line.find('=')
                new_line = 'name' + str(i) + line[eq_idx:]
                i += 1
                f.write(new_line)
            else:
                f.write(line)
    return common_modification_list


# 读取盲搜结果并生成未知质量数修饰
def mass_diff_list_generate(res_path, current_path):
    summary_path = os.path.join(res_path, 'pFind.summary')
    with open(summary_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    mass_diff_lines = []
    i = 0
    while i < len(lines):
        if 'Modifications:' in lines[i]:
            i += 1
            mass_diff_lines.append(lines[i])
            i += 1
            j = 0
            while j < 20 and i < len(lines):
                if 'PFIND' in lines[i]:
                    mass_diff_lines.append(lines[i])
                    j += 1
                i += 1
        i += 1
    # 写入txt文件方便筛选
    mass_diff_path = os.path.join(current_path, 'mass_diff_list.txt')
    with open(mass_diff_path, 'w', encoding='utf-8') as f:
        for line in mass_diff_lines:
            f.write(line)


