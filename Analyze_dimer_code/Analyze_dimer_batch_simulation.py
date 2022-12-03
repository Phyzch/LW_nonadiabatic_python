from matplotlib import pyplot, gridspec
from util import *
from Analyze_dimer_code.Analyze_dimer_batch_siimulation_auxiliary_func import *


def plot_dimer_survival_prob_selected_states_batch_simulation():
    '''

    :param file_path:
    :return:
    '''
    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/output_file/Vt=50/"
    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, survival_prob_list, time_list = Read_dimer_survival_prob_all_states(file_path)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    selected_index = 0
    label = "|v1> = " + str(monomer1_quantum_num[selected_index]) + " |v2> = " + str(monomer2_quantum_num[selected_index])

    ax.plot( time_list, survival_prob_list[selected_index] , linewidth = 2, label = label )

    ax.legend(loc = 'best')
    ax.set_yscale('log')

    plt.show()

def analyze_dilution_factor_and_T_T_prime_phase_diagram():
    '''

    :return:
    '''
    # anharmonic coupling
    V0 = 300
    scaling_factor = 0.3

    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/output_file/"

    file_path1 = "Vt=0"
    file_path2 = "Vt=50"
    file_path3 = "Vt=100"
    file_path4 = "Vt=200"
    file_path5 = "Vt=363"
    file_path6 = "Vt=500"

    ## additional results for states with low energy
    file_path1_1 = "Vt=0_low_energy"
    file_path2_1 = "Vt=50_low_energy"
    file_path3_1 = "Vt=100_low_energy"
    file_path4_1 = "Vt=200_low_energy"
    file_path5_1 = "Vt=363_low_energy"
    file_path6_1 = "Vt=500_low_energy"

    file_path_list = [file_path1, file_path2, file_path3, file_path4, file_path5, file_path6,
                      file_path1_1, file_path2_1, file_path3_1, file_path4_1, file_path5_1, file_path6_1]

    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    Vt_list = [0, 50, 100, 200, 363, 500, 0, 50, 100, 200, 363, 500]
    frequency_list_monomer1 = np.array([890, 727, 345, 1117, 1158])
    frequency_list_monomer2 = np.array([890, 727, 345, 1117, 1158])

    monomer1_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])
    monomer2_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])

    dimer_alpha = [monomer1_alpha , monomer2_alpha]
    dimer_frequency = [ frequency_list_monomer1 , frequency_list_monomer2 ]

    analyze_dilution_factor_and_T_T_prime_phase_diagram_subroutine(folder_path, file_path_list, Vt_list,
                                                                   V0, scaling_factor, dimer_alpha,
                                                                   dimer_frequency, save_bool)

def analyze_dilution_factor_and_T_T_prime_phase_diagram_subroutine(folder_path, file_path_list, Vt_list, V0, scaling_factor, dimer_alpha, dimer_frequency, save_bool):
    '''

    :param folder_path:
    :param file_path_list:
    :param Vt_list:
    :param V0:
    :param scaling_factor:
    :param dimer_alpha:
    :param dimer_frequency:
    :param save_bool:
    :return:
    '''

    monomer1_alpha = dimer_alpha[0]
    monomer2_alpha = dimer_alpha[1]

    monomer1_frequency = dimer_frequency[0]
    monomer2_frequency = dimer_frequency[1]

    path_num = len(file_path_list)

    T_list = []
    T_prime_list = []
    dilution_factor_list = []

    for i in range(path_num):
        file_path = file_path_list[i]
        Vt = Vt_list[i]
        state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor, time_list = compute_dimer_state_dilution_factor(file_path)

        dilution_factor = dilution_factor.tolist()

        # transition factor in same PES
        T = estimate_transition_factor_T_dimer_lists(V0, scaling_factor, file_path)
        T = T.tolist()

        # transition factor in different PES
        T_prime = estimate_nonadiabatic_transition_factor_T_prime_dimer_list(dimer_alpha, dimer_frequency, file_path, Vt)
        T_prime = T_prime.tolist()

        # concatenate list
        T_list = T_list + T
        T_prime_list = T_prime_list + T_prime
        dilution_factor_list = dilution_factor_list + dilution_factor

    dilution_factor_list = np.array(dilution_factor_list)
    T_list = np.array(T_list)
    T_prime_list = np.array(T_prime_list)

    dilution_factor_criterion1 = 0.1
    dilution_factor_criterion2 = 0.01
    length = len(dilution_factor_list)

    localization_index_list = [index for index in range(length) if
                               dilution_factor_list[index] > dilution_factor_criterion1]
    mixing_index_list = [index for index in range(length) if
                         dilution_factor_list[index] <= dilution_factor_criterion1 and dilution_factor_list[
                             index] > dilution_factor_criterion2]

    extended_state_index_list = [index for index in range(length) if
                                 dilution_factor_list[index] <= dilution_factor_criterion2]

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    ax.scatter(T_list[localization_index_list] , T_prime_list[localization_index_list] , marker = 'o' , color = 'red' , s=40)
    ax.scatter(T_list[extended_state_index_list] , T_prime_list[extended_state_index_list] , marker =  'v', color = 'blue', s=40)
    ax.scatter(T_list[mixing_index_list], T_prime_list[mixing_index_list], marker = 's' , color = 'purple', s = 40 )

    T_boundary = np.linspace(0,1, 100)
    T_prime_boundary = 1 - T_boundary
    ax.plot(T_boundary, T_prime_boundary, color = 'black')

    ax.set_xlabel("$T$")
    ax.set_ylabel(" $T'$ ")
    ax.set_ylim([-0.05, 0.6 ])
    ax.set_xlim([0,1])

    if save_bool:
        fig_name = "dilution factor vs T.svg"
        fig_name = os.path.join(folder_path, fig_name)
        fig.savefig(fig_name)

    plt.show()