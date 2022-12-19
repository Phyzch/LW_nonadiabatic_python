from util import *
from Analyze_dimer_code.Analyze_dimer_batch_siimulation_survival_prob_auxiliary_func import *
from estimate_transition_factor_T_full_mode import  estimate_transition_energy
from scipy.optimize import root
from estimate_transition_factor_T_full_mode import estimate_transition_factor_with_freq_and_energy_dimer
from estimate_transition_factor_T_prime_full_mode import estimate_T_prime_prefactor_subroutine_dimer

def analyze_dilution_factor_and_E_Vt_phase_diagram_random_Vt():
    '''

    :return:
    '''
    V0 = 3050
    frequency_list =  np.array([890, 727, 345, 1117, 1158])
    dof = len(frequency_list)

    scaling_factor_list = np.sqrt(frequency_list) / 270
    scaling_factor = np.prod( np.power(scaling_factor_list , 1/dof) )

    # estimate using equal energy per monomer assumption
    alpha_list = np.array([0.169, 0.163, 0.127, 0.101, 0.101])

    save_bool = True
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/" \
                  "BChl_dimer_model/5_mode/batch_simulation_Bigwood_scaling/batch_simulation_random_Vt/output/"

    file_path1 = "Vt=0-100"
    file_path2 = "Vt=100-200"
    file_path3 = "Vt=200-300"
    file_path4 = "Vt=300-400"


    file_path_list = [file_path1, file_path2, file_path3, file_path4]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    color_list = ['blue', 'orange', 'gray', 'brown', 'black'  ]


    analyze_dilution_factor_and_E_Vt_phase_diagram_random_Vt_subroutine(folder_path, file_path_list, frequency_list, scaling_factor, alpha_list, save_bool)

def Read_Vt_coupling_info_dimer(folder_path):
    '''

    :param file_path:
    :return:
    '''
    file_name = "sampling_state_Vt.txt"
    file_path = os.path.join(folder_path, file_name)

    monomer1_number_list = []
    monomer2_number_list = []

    with open (file_path) as f:
        data = f.read()
        data = re.split('\n', data)
        data = [i for i in data if i != '']
        datalen = len(data)

        line = data[0]
        line = re.split(' ', line)
        state_num = int(line[0])

        line = data[1]
        line = re.split(' ' , line)
        state_energy = [float(i) for i in line if i!= '']

        line = data[2]
        line = re.split(' ' , line)
        Vt_list = [float(i) for i in line if i!= '']

        line_index = 3

        for i in range(state_num):
            line = data[line_index]
            line = re.split(' ', line)
            mode_number = [int(i) for i in line if i!= '']
            nmode = int( len(mode_number) / 2 )

            monomer1 = mode_number[:nmode]
            monomer2 = mode_number[nmode:]

            monomer1_number_list.append(monomer1)
            monomer2_number_list.append(monomer2)

            line_index = line_index + 1

    state_energy = np.array(state_energy)
    Vt_list = np.array(Vt_list)

    monomer1_number_list = np.array(monomer1_number_list)
    monomer2_number_list = np.array(monomer2_number_list)

    return state_num, state_energy, Vt_list, monomer1_number_list, monomer2_number_list

def analyze_dilution_factor_and_E_Vt_phase_diagram_random_Vt_subroutine(fig_path, file_path_list, monomer_frequency, scaling_factor, alpha_list,  save_bool):
    '''

    :param fig_path:
    :param file_path_list:
    :param dimer_alpha:
    :param dimer_frequency:
    :param save_bool:
    :return:
    '''

    path_num = len(file_path_list)

    state_energy_list = []
    Vt_for_state_list = []
    dilution_factor_list = []

    for i in range(path_num):
        file_path = file_path_list[i]
        state_num, state_energy, Vt_list_in_one_file, monomer1_number_list, monomer2_number_list = Read_Vt_coupling_info_dimer(file_path)
        Vt_for_state_list = Vt_for_state_list + Vt_list_in_one_file.tolist()

        state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor, time_list = compute_dimer_state_dilution_factor(file_path)
        state_energy_list = state_energy_list + state_energy.tolist()
        dilution_factor_list = dilution_factor_list + dilution_factor.tolist()

    dilution_factor_list = np.array(dilution_factor_list)
    state_energy_list = np.array(state_energy_list)
    Vt_for_state_list = np.array(Vt_for_state_list)

    Vt_cutoff = 100
    small_Vt_state_index = [index for index in range(len(Vt_for_state_list)) if Vt_for_state_list[index] < Vt_cutoff]
    dilution_factor_small_Vt = dilution_factor_list[small_Vt_state_index]
    dilution_factor_criterion = np.sqrt(np.min(dilution_factor_small_Vt))

    # dilution_factor_criterion = 0.04
    print("dilution factor criterion: " + str(dilution_factor_criterion))

    length = len(dilution_factor_list)
    localization_index_list = [index for index in range(length) if
                               dilution_factor_list[index] > dilution_factor_criterion]

    extended_state_index_list = [index for index in range(length) if
                                 dilution_factor_list[index] <= dilution_factor_criterion]

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    ax.scatter(state_energy_list[localization_index_list], Vt_for_state_list[localization_index_list], marker='o',
               color='red', s=40)
    ax.scatter(state_energy_list[extended_state_index_list], Vt_for_state_list[extended_state_index_list], marker='v',
               color='black', s=40)
    # ax.scatter(T_list[mixing_index_list], T_prime_list[mixing_index_list], marker = 's' , color = 'purple', s = 40 )

    ax.set_xlabel("$E$")
    ax.set_ylabel(" $V_{t}$ ")

    # draw line on E, Vt phase diagram about T_FR + T_NA = 1:
    result = root(estimate_transition_energy, x0 = np.array([15000]) , args=(monomer_frequency,0, scaling_factor))
    critical_energy = result.x[0]
    min_energy_sample = 10000
    energy_data_num = 50
    energy_list = np.linspace(min_energy_sample,critical_energy, energy_data_num)
    # equal energy assumption
    critical_Vt_list = []
    for i in range(energy_data_num):
        energy = energy_list[i]
        T_prime_prefactor = estimate_T_prime_prefactor_subroutine_dimer(alpha_list, monomer_frequency, energy, equal_energy_assumption= True)
        T = estimate_transition_factor_with_freq_and_energy_dimer(energy, monomer_frequency,scaling_factor, equal_energy_assumption= True)
        critical_Vt = (1 - T)/T_prime_prefactor
        critical_Vt_list.append(critical_Vt)

    ax.plot(energy_list, critical_Vt_list, color = 'black' , linewidth = 3)
    ax.set_ylim([0,410])
    ax.set_xlim([np.min(state_energy_list) , 25000])

    if save_bool:
        fig_name = "dilution factor vs energy and Vt.svg"
        fig_name = os.path.join(fig_path, fig_name)
        fig.savefig(fig_name)

    plt.show()