from util import *
import random

def generate_sampling_state_info_file():
    '''

    :param folder_path:
    :return:
    '''
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/try/"
    file_name = "sampling_state_info.txt"
    file_path = os.path.join(folder_path, file_name)

    sampling_state_num = 10
    monomer_mode_num = 5

    nmax_sampling_list = np.array([ 3, 3, 3, 3 , 3 ])
    upper_energy_cutoff = 5000
    lower_energy_cutoff = 3000

    frequency_list = np.array([890, 727, 345, 1117, 1158])

    dimer_state_list = generate_random_sampling_state(sampling_state_num, monomer_mode_num, nmax_sampling_list, frequency_list, lower_energy_cutoff, upper_energy_cutoff)

    output_sampling_state_info_to_file(file_path, sampling_state_num, monomer_mode_num, dimer_state_list )


def generate_random_sampling_state(sampling_state_num, monomer_mode_num, nmax_sampling_list, frequency_list , lower_energy_cutoff, upper_energy_cutoff):
    '''
    energy of generated random state should be smaller than energy_cutoff
    nmax_list : nmax for each dof. (maximum quantum num)
    :param sampling_state_num:
    :param monomer_mode_num:
    :param nmax_list:
    :param energy_cutoff:
    :return:
    '''

    dimer_state_list = []

    while len(dimer_state_list) < sampling_state_num:
        dimer_state = []
        monomer_state1 = []
        for j in range(monomer_mode_num):
            monomer_state1.append( round( nmax_sampling_list[j] * random.random() ) )

        monomer_state2 = []
        for j in range(monomer_mode_num):
            monomer_state2.append( round( nmax_sampling_list[j] * random.random()  ) )

        dimer_state = [monomer_state1, monomer_state2]

        energy = 0
        energy = energy + np.sum( np.array(monomer_state1) * frequency_list )
        energy = energy + np.sum( np.array(monomer_state2) * frequency_list )

        if energy > upper_energy_cutoff or energy < lower_energy_cutoff:
            continue

        exist = search_list(dimer_state_list, dimer_state)
        if exist:
            continue

        dimer_state_list.append(dimer_state)

    return dimer_state_list

def output_sampling_state_info_to_file(file_path, sampling_state_num, monomer_mode_num, dimer_state_list ):
    '''

    :param file_path:
    :param sampling_state_num:
    :param monomer_mode_num:
    :param dimer_state_list:
    :return:
    '''
    with open(file_path, "w") as f:
        f.write(str(sampling_state_num) + " ")
        f.write(str(monomer_mode_num) + " ")
        f.write("\n")

        for state_index in range(sampling_state_num):
            monomer1 = dimer_state_list[state_index][0]
            monomer2 = dimer_state_list[state_index][1]
            for mode_index in range(monomer_mode_num):
                f.write(str(monomer1[mode_index]) + " ")
            f.write("     ")
            for mode_index in range(monomer_mode_num):
                f.write(str(monomer2[mode_index]) + " ")
            f.write("\n")

def search_list (list, element):
    for i in range(len(list)):
        if(element == list[i]):
            return True

    return False


