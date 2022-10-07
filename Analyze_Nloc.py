'''
analyze N_{loc} for simulation.  results store in diff_PES_Nloc.txt and same_PES_Nloc.txt
'''
from util import *

def read_Nloc(file_name):
    '''

    :param file_name:
    :return:
    '''
    with open(file_name) as f:
        data = f.read().splitlines()
        data = [i for i in data if i != '']
        datalen = len(data)

        line = re.split(' ', data[0])
        line = [int(i) for i in line]
        state_num , mode_num = line

        mode_index_list = []
        Nloc_list = []

        line_index = 0
        for i in range(int(state_num)):
            line_index = line_index + 1
            line = re.split(' ', data[line_index])
            mode_index = [int(i) for i in line if i!= '']

            line_index = line_index + 1
            Nloc = float(data[line_index])

            mode_index_list.append(mode_index)
            Nloc_list.append(Nloc)

        mode_index_list = np.array(mode_index_list)
        Nloc_list = np.array(Nloc_list)

        return state_num , mode_num, mode_index_list, Nloc_list

def read_Nloc_anharmonic_nonadiabatic_coupling(folder_path):
    '''

    :param folder_path:
    :return:
    '''
    # coupling within same PES
    same_PES_file_name = "same_PES_Nloc.txt"
    same_PES_file_path = os.path.join(folder_path, same_PES_file_name)

    diff_PES_file_name = "diff_PES_Nloc.txt"
    diff_PES_file_path = os.path.join(folder_path, diff_PES_file_name)

    state_num, mode_num, mode_index_list, same_PES_Nloc_list = read_Nloc(same_PES_file_path)

    _, _, _, diff_PES_Nloc_list = read_Nloc(diff_PES_file_path)

    return state_num, mode_num, mode_index_list, same_PES_Nloc_list, diff_PES_Nloc_list