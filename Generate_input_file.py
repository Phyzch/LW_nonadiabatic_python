import numpy as np
import os
import re
from shutil import copy

def Generate_input_files(folder_path):
    sampling_state_index_file_path = os.path.join(folder_path, "sampling_state_info.txt")
    input_file_path = os.path.join(folder_path, "input.txt")

    tunneling_strength = [0, 5, 10 ,20, 30, 50, 70, 100 , 200 , 300 ]
    scaling_factor = [ 0.16, 0.17 ]

    tunneling_num = len(tunneling_strength)
    scaling_num = len(scaling_factor)

    if (os.path.exists(sampling_state_index_file_path) and os.path.exists(input_file_path)):
        for i in range(tunneling_num):
            subfolder_path = os.path.join(folder_path, str(tunneling_strength[i]))
            if(not os.path.exists(subfolder_path)):
                os.mkdir(subfolder_path)
            for j in range(scaling_num):
                subsubfolder_path = os.path.join(subfolder_path, str(scaling_factor[j]) )
                if(not os.path.exists(subsubfolder_path)):
                    os.mkdir(subsubfolder_path)

                # copy files
                input_subfolder_path_copy = os.path.join(subsubfolder_path,"input(copy).txt")
                input_subfolder_path = os.path.join(subsubfolder_path,"input.txt")
                sampling_state_index_subfolder_file_path = os.path.join(subsubfolder_path, "sampling_state_info.txt")

                copy(input_file_path, input_subfolder_path_copy)
                copy(sampling_state_index_file_path, sampling_state_index_subfolder_file_path)

                # change input file content
                Change_input_file(input_subfolder_path_copy, input_subfolder_path, tunneling_strength[i], scaling_factor[j])

                os.remove(input_subfolder_path_copy)


    else:
        raise NameError("sampling_state_info.txt or input.txt does not exist.")

def Change_input_file(input_path, output_path, tunneling_strength, scaling_factor):
    fin = open(input_path, "rt")
    fout = open(output_path, "wt")

    data= fin.readlines()

    # read and write first line
    fout.write(data[0])

    line = data[1]
    line = re.split(' ', line)
    line = [float(i) for i in line if i!='']

    # replace scaling factor and tunneling_strength to value we want
    line[7] = scaling_factor
    line[12] = tunneling_strength
    # Rmax should be int
    line[5] = int(line[5])
    for i in range(len(line)):
        fout.write(str(line[i]) + "  ")

    fout.write("\n")

    line_index = 2
    datalen = len(data)
    for i in range(line_index, datalen):
        fout.write(data[i])

    fin.close()
    fout.close()