# 对盲搜鉴定到的修饰，进行质量校正
import os 
import numpy as np 
from collections import Counter 


element_dict={
    "C": 12.0000000,
    "H": 1.0078246,
    "Pm": 1.00727647012,
    "N": 14.0030732,
    "O": 15.9949141,
    "S": 31.972070
}

amino_acid_dict={
    "A" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "C" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*1,
    "D" : element_dict["C"]*4 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*3 + element_dict["S"]*0,
    "E" : element_dict["C"]*5 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*3 + element_dict["S"]*0,
    "F" : element_dict["C"]*9 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "G" : element_dict["C"]*2 + element_dict["H"]*3 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "H" : element_dict["C"]*6 + element_dict["H"]*7 + element_dict["N"]*3 + element_dict["O"]*1 + element_dict["S"]*0,
    "I" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "K" : element_dict["C"]*6 + element_dict["H"]*12 + element_dict["N"]*2 + element_dict["O"]*1 + element_dict["S"]*0,
    "L" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "M" : element_dict["C"]*5 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*1,
    "N" : element_dict["C"]*4 + element_dict["H"]*6 + element_dict["N"]*2 + element_dict["O"]*2 + element_dict["S"]*0,
    "P" : element_dict["C"]*5 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "Q" : element_dict["C"]*5 + element_dict["H"]*8 + element_dict["N"]*2 + element_dict["O"]*2 + element_dict["S"]*0,
    "R" : element_dict["C"]*6 + element_dict["H"]*12 + element_dict["N"]*4 + element_dict["O"]*1 + element_dict["S"]*0,
    "S" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
    "T" : element_dict["C"]*4 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
    "V" : element_dict["C"]*5 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "W" : element_dict["C"]*11 + element_dict["H"]*10 + element_dict["N"]*2 + element_dict["O"]*1 + element_dict["S"]*0,
    "X" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "Y" : element_dict["C"]*9 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
}

h2o_mass = element_dict["H"]*2 + element_dict["O"]*1
proton_mass = element_dict['Pm']


# 读取选择的常见修饰，返回dict
def common_dict_create(current_path):
    modification_path = os.path.join(current_path, 'modification-null.ini')
    with open(modification_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    i = 1
    common_dict = {} 
    while i < len(lines):
        if len(lines[i]) < 2:
            break
        mod_name = lines[i].split()[0]
        eq_idx = mod_name.find('=')
        mod_name = mod_name[eq_idx+1:]
        mod_mass = lines[i+1].split()[2]
        common_dict[mod_name] = float(mod_mass)
        i += 2
    return common_dict


# 计算系统误差
def system_shift_compute(lines):
    mass_shift = []
    for line in lines:
        line = line.split('\t')
        if len(line[10]):
            continue
        else:
            mass_shift.append(float(line[7]))
    system_shift = np.mean(mass_shift)
    return system_shift


# 计算修饰的精确质量 
def accurate_mass_compute(lines, mass, common_dict):
    mass_list = []
    for line in lines:
        if mass not in line:
            continue
        line = line.split('\t')
        mod_list = line[10].split(';')[:-1]
        parent_mass = float(line[2])
        sequence = line[5]
        amino_mass = 0.0
        for a in sequence:
            amino_mass += amino_acid_dict[a]
        mod_mass = parent_mass - amino_mass - proton_mass - h2o_mass
        if len(mod_list) > 1:
            for mod in mod_list:
                if mass in mod:
                    continue
                mod = mod.split(',')[1]
                mod_mass -= common_dict[mod]
        mass_list.append(mod_mass)
    return np.mean(mass_list)



def mass_correct(current_path, blind_path, mass_diff_list):
    # 读取常见修饰列表
    common_dict = common_dict_create(current_path)

    # 读取盲搜鉴定结果文件
    with open(blind_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    lines = lines[1:]

    # 计算系统误差
    system_shift = system_shift_compute(lines)
    
    # 计算未知修饰的精确质量
    mass_dict = {}
    for mass_diff in mass_diff_list:
        accurate_mass_diff = accurate_mass_compute(lines, mass_diff, common_dict)
        mass_dict[mass_diff] = float('%.6f'%(accurate_mass_diff - system_shift))
    return mass_dict



# 删除质量小于200Da的偏差和负数的修饰 
def small_delta_filter(mass_difference_list, min_mass_modification):
    name2mass = {}
    new_mass_diff_list = []
    for mass_diff in mass_difference_list:
        mass = mass_diff.split('_')[-1]
        if mass[0] == '-':
            continue
        mass = float(mass)
        if mass < min_mass_modification:
            continue
        name2mass[mass_diff] = mass
        new_mass_diff_list.append(mass_diff)
    return name2mass, new_mass_diff_list  


# 筛选修饰质量差满足设定的插值的修饰 
# mass_diff_diff = 6.0201 
# 应该是20ppm，但是实际太小完全没有删选度
def mass_diff_diff_filter(name2mass, mass_diff_list, mass_diff_diff):
    refined_list = []
    for i in range(len(mass_diff_list)):
        for j in range(i+1, len(mass_diff_list)):
            mass_left = name2mass[mass_diff_list[i]]
            mass_right = name2mass[mass_diff_list[j]]
            delta_mass = abs(mass_right -  mass_left)
            if abs(delta_mass - mass_diff_diff) < 0.1:
                refined_list.append(mass_diff_list[i])
                refined_list.append(mass_diff_list[j])
    new_refined_list = []
    for mod_name in refined_list:
        flag_rep = False 
        for ref_mod in new_refined_list:
            if abs(name2mass[mod_name] - name2mass[ref_mod]) < 0.05:
                flag_rep = True
        if flag_rep == False:
            new_refined_list.append(mod_name)
    return new_refined_list 


# 统计修饰发生的位点分布 
def mass_static(blind_path, current_path, mass_diff_list): 
    mod_position_dict = {} 
    mod_number_dict = {} 
    for mass in mass_diff_list:
        mod_position_dict[mass] = [] 
        mod_number_dict[mass] = 0 
    
    with open(blind_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    spectra_num = len(lines)
    for i in range(1, spectra_num):
        if len(lines[i]) < 4:
            break
        sequence = lines[i].split('\t')[5]
        mod_list = lines[i].split('\t')[10].split(';')[:-1]
        for mod in mod_list:
            pos, mod_name = mod.split(',')
            if mod_name in mass_diff_list:
                pos = int(pos) 
                mod_number_dict[mod_name] += 1
                if pos == 0 or pos == 1:
                    mod_position_dict[mod_name].append(sequence[0])
                    mod_position_dict[mod_name].append('N-SIDE')
                elif pos >= len(sequence):
                    mod_position_dict[mod_name].append(sequence[-1])
                    mod_position_dict[mod_name].append('C_SIDE')
                else:
                    mod_position_dict[mod_name].append(sequence[pos-1])
    
    mod_static_dict = {}
    for mod_name in mass_diff_list:
        counter = Counter(mod_position_dict[mod_name])
        mod_static_dict[mod_name] = counter 
    
    return mod_static_dict, mod_number_dict  


# 写pchem.summary结果文件
def summary_write(current_path, mod_static_dict, mod_number_dict): 
    lines = []
    lines.append('rank \tmass modification \tTotal number \tposition: occurance times\n')
    i = 1 
    for mod in mod_static_dict:
        line = str(i) + '\t' + mod + '\t' + str(mod_number_dict[mod]) + '\t' + str(mod_static_dict[mod])[8:-1] + '\n'
        lines.append(line)
        i += 1
    with open(os.path.join(current_path, 'pChem.summary'), 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line)
    # print(mod_static_dict)


# 选择不重复的topk来做限定式搜索
def mass_select(mass_diff_list, k, name2mass):
    selected_list = []
    i = 0 
    for mod in mass_diff_list: 
        flag = True
        for reference in selected_list:
            if abs(name2mass[mod]-name2mass[reference]) < 0.03:
                flag = False
        if flag == True:
            selected_list.append(mod)
            i += 1 
        if i >= k:
            break
    return selected_list 
        

