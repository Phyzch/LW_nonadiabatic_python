from Analyze_dimer_code.Analyze_dimer_batch_siimulation_survival_prob_auxiliary_func import *


def plot_dimer_survival_prob_selected_states_batch_simulation():
    '''

    :param file_path:
    :return:
    '''
    save_bool = False
    fig_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/test_anharmonic_effect/"

    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/test_anharmonic_effect/V3=0"
    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, survival_prob_list, time_list = Read_dimer_survival_prob_all_states(file_path)

    file_path2 = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/test_anharmonic_effect/V3=100,a=0.25"
    state_num, nmode, state_energy2, monomer1_quantum_num2, monomer2_quantum_num2, survival_prob_list2, time_list2 = Read_dimer_survival_prob_all_states(file_path2)

    file_path3 = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/test_anharmonic_effect/V3=300,a=0.25"
    state_num, nmode, state_energy3, monomer1_quantum_num3, monomer2_quantum_num3, survival_prob_list3, time_list3 = Read_dimer_survival_prob_all_states(file_path3)


    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    selected_index = 0
    label = "|v1> = " + str(monomer1_quantum_num[selected_index]) + " |v2> = " + str(monomer2_quantum_num[selected_index]) + " a=0"
    label2 = "|v1> = " + str(monomer1_quantum_num2[selected_index]) + " |v2> = " + str(
        monomer2_quantum_num2[selected_index]) + " V3=100 a=0.25"
    label3 = "|v1> = " + str(monomer1_quantum_num2[selected_index]) + " |v2> = " + str(
        monomer2_quantum_num2[selected_index]) + " V3=300, a=0.25"

    final_time = 2
    dt = time_list[1] - time_list[0]
    final_time_index = int(final_time / dt)

    ax.plot( time_list[:final_time_index], survival_prob_list[selected_index, : final_time_index] , linewidth = 2, label = label )
    ax.plot(time_list2[:final_time_index], survival_prob_list2[selected_index, : final_time_index], linewidth = 2, label = label2)
    ax.plot(time_list3[:final_time_index], survival_prob_list3[selected_index, : final_time_index], linewidth = 2, label = label3)


    ax.legend(loc = 'best')
    ax.set_yscale('log')

    if save_bool:
        fig_name = "survival prob compare.svg"
        fig_path = os.path.join(fig_path, fig_name)
        fig.savefig(fig_path)

    plt.show()

def analyze_dilution_factor_and_T_T_prime_phase_diagram_main():
    '''

    :return:
    '''
    # V0 = 300, a = 0.3
    # analyze_dilution_factor_and_T_T_prime_phase_diagram_V0_300()

    # V0 = 3050,  a = geometric mean of 'a' for each mode
    analyze_dilution_factor_and_T_T_prime_phase_diagram_V0_3050()

def analyze_dilution_factor_and_T_T_prime_phase_diagram_V0_3050():
    '''
    Bigwood 1998 PNAS formula
    :return:
    '''
    V0 = 3050
    frequency_list =  np.array([890, 727, 345, 1117, 1158])
    dof = len(frequency_list)

    scaling_factor_list = np.sqrt(frequency_list) / 270
    scaling_factor = np.prod( np.power(scaling_factor_list , 1/dof) )

    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation_Bigwood_scaling/batch_simulation_output/"

    file_path1 = "Vt=0"
    file_path2 = "Vt=50"
    file_path3 = "Vt=100"
    file_path4 = "Vt=200"
    file_path5 = "Vt=363"


    file_path_list = [file_path1, file_path2, file_path3, file_path4, file_path5]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    Vt_list = [0, 50, 100, 200, 363]
    color_list = ['blue', 'orange', 'gray', 'brown', 'black'  ]

    frequency_list_monomer1 = np.array([890, 727, 345, 1117, 1158])
    frequency_list_monomer2 = np.array([890, 727, 345, 1117, 1158])

    monomer1_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])
    monomer2_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])

    dimer_alpha = [monomer1_alpha , monomer2_alpha]
    dimer_frequency = [ frequency_list_monomer1 , frequency_list_monomer2 ]

    label_list = ['Vt=0' , 'Vt=50', 'Vt=100', 'Vt=200', 'Vt=363']

    analyze_dilution_factor_and_T_T_prime_log_scale_plot(folder_path, file_path_list, Vt_list,
                                                         V0, scaling_factor, dimer_alpha, dimer_frequency,
                                                         save_bool, color_list, label_list)

    # analyze_dilution_factor_and_T_T_prime_phase_diagram_subroutine(folder_path, file_path_list, Vt_list,
    #                                                                V0, scaling_factor, dimer_alpha, dimer_frequency, save_bool)

    # analyze_dilution_factor_and_E_Vt_prime_phase_diagram_subroutine(folder_path, file_path_list,
    #                                                                 Vt_list, dimer_alpha, dimer_frequency, save_bool)

def analyze_dilution_factor_and_T_T_prime_phase_diagram_V0_300():
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

    ## additioinal_results_for_states_with_middle_energy
    file_path1_2 = "Vt=0_middle_energy"
    file_path2_2 = "Vt=50_middle_energy"
    file_path3_2 = "Vt=100_middle_energy"
    file_path4_2 = "Vt=200_middle_energy"
    file_path5_2 = "Vt=363_middle_energy"
    file_path6_2 = "Vt=500_middle_energy"

    ## additioinal_results_for_states_with_middle_energy
    file_path1_3 = "Vt=0_high_energy"
    file_path2_3 = "Vt=50_high_energy"
    file_path3_3 = "Vt=100_high_energy"
    file_path4_3 = "Vt=200_high_energy"
    file_path5_3 = "Vt=363_high_energy"
    file_path6_3 = "Vt=500_high_energy"

    file_path_list = [file_path1,  file_path3, file_path4, file_path5, file_path6,
                      file_path1_1,  file_path3_1, file_path4_1, file_path5_1, file_path6_1,
                      file_path1_2,  file_path3_2, file_path4_2, file_path5_2, file_path6_2,
                      file_path1_3, file_path3_3, file_path4_3, file_path5_3, file_path6_3]

    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    Vt_list = [0,  100, 200, 363, 500] * 4
    color_list = ['blue', 'orange', 'red', 'brown', 'black' ] * 4

    frequency_list_monomer1 = np.array([890, 727, 345, 1117, 1158])
    frequency_list_monomer2 = np.array([890, 727, 345, 1117, 1158])

    monomer1_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])
    monomer2_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])

    dimer_alpha = [monomer1_alpha , monomer2_alpha]
    dimer_frequency = [ frequency_list_monomer1 , frequency_list_monomer2 ]

    # analyze_dilution_factor_and_T_T_prime_phase_diagram_subroutine(folder_path, file_path_list, Vt_list,
    #                                                                V0, scaling_factor, dimer_alpha,
    #                                                                dimer_frequency, save_bool)

    analyze_dilution_factor_and_T_T_prime_log_scale_plot(folder_path, file_path_list, Vt_list, V0, scaling_factor,
                                                         dimer_alpha, dimer_frequency, save_bool , color_list)



def analyze_dilution_factor_and_T_T_prime_log_scale_plot(folder_path, file_path_list, Vt_list, V0, scaling_factor, dimer_alpha, dimer_frequency, save_bool, color_list, label_list):
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

    T_T_prime_sum_list = []
    dilution_factor_list = []

    for i in range(path_num):
        file_path = file_path_list[i]
        Vt = Vt_list[i]
        state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor, time_list = compute_dimer_state_dilution_factor(
            file_path)

        dilution_factor = dilution_factor.tolist()

        # transition factor in same PES
        T = estimate_transition_factor_T_dimer_lists(V0, scaling_factor, file_path)

        # transition factor in different PES
        T_prime = estimate_nonadiabatic_transition_factor_T_prime_dimer_list(dimer_alpha, dimer_frequency, file_path,
                                                                             Vt)

        T_T_prime_sum = np.array(T) + np.array(T_prime)

        # concatenate list
        T_T_prime_sum_list.append(T_T_prime_sum)
        dilution_factor_list.append(dilution_factor)

    dilution_factor_list = np.array(dilution_factor_list)
    T_T_prime_sum_list = np.array(T_T_prime_sum_list)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    T_T_prime_sum_list_len = len(T_T_prime_sum_list)
    for i in range(T_T_prime_sum_list_len):
        ax.scatter(T_T_prime_sum_list[i], dilution_factor_list[i], marker = 'o' , color = color_list[i] , label = label_list[i], s=30 )

    ax.set_xlabel(" T + T' ")
    ax.set_ylabel('$\sigma$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc = 'best')
    ax.set_xlim([0.1, np.max(T_T_prime_sum_list) + 0.5 ])

    if save_bool:
        fig_name = "dilution factor vs T+T' log scale.svg"
        fig_name = os.path.join(folder_path, fig_name)
        fig.savefig(fig_name)

    # for each Vt plot figure

    # for i in range(path_num):
    #     figi = plt.figure(figsize=(10, 10))
    #     speci = gridspec.GridSpec(nrows=1, ncols=1, figure=figi)
    #     axi =  figi.add_subplot(speci[0,0])
    #
    #     axi.scatter(T_T_prime_sum_list[i], dilution_factor_list[i], marker='o', color=color_list[i], label=label_list[i],
    #                s=30)
    #     axi.set_xscale('log')
    #     axi.set_yscale('log')
    #     axi.set_xlabel(" T + T' ")
    #     axi.set_ylabel('$\sigma$')
    #     axi.legend(loc = 'best')
    #     axi.set_xlim([0.1, np.max(T_T_prime_sum_list) + 0.5 ])
    #     axi.set_ylim([ 5* pow(10,-5) , 1 ])
    #
    #     if save_bool:
    #         fig_name = "dilution factor vs T+T' log scale" + str(i) +".svg"
    #         fig_name = os.path.join(folder_path, fig_name)
    #         figi.savefig(fig_name)

    plt.show()


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


    dilution_factor_criterion2 = 0.01
    length = len(dilution_factor_list)

    localization_index_list = [index for index in range(length) if
                               dilution_factor_list[index] > dilution_factor_criterion2]

    # mixing_index_list = [index for index in range(length) if
    #                      dilution_factor_list[index] <= dilution_factor_criterion1 and dilution_factor_list[
    #                          index] > dilution_factor_criterion2]

    extended_state_index_list = [index for index in range(length) if
                                 dilution_factor_list[index] <= dilution_factor_criterion2]

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    ax.scatter(T_list[localization_index_list] , T_prime_list[localization_index_list] , marker = 'o' , color = 'red' , s=40)
    ax.scatter(T_list[extended_state_index_list] , T_prime_list[extended_state_index_list] , marker =  'v', color = 'black', s=40)
    # ax.scatter(T_list[mixing_index_list], T_prime_list[mixing_index_list], marker = 's' , color = 'purple', s = 40 )

    T_boundary = np.linspace(0,1, 100)
    T_prime_boundary = 1 - T_boundary
    ax.plot(T_boundary, T_prime_boundary, color = 'black' , linewidth = 3 )

    ax.set_xlabel("$T$")
    ax.set_ylabel(" $T'$ ")
    ax.set_ylim([-0.05, 0.7 ])
    ax.set_xlim([0, 3.8 ])

    if save_bool:
        fig_name = "dilution factor vs T.svg"
        fig_name = os.path.join(folder_path, fig_name)
        fig.savefig(fig_name)

    plt.show()


def analyze_dilution_factor_and_E_Vt_prime_phase_diagram_subroutine(folder_path, file_path_list, Vt_list, dimer_alpha, dimer_frequency, save_bool):
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

    state_energy_list = []
    Vt_for_state_list = []
    dilution_factor_list = []

    for i in range(path_num):
        file_path = file_path_list[i]
        Vt = Vt_list[i]
        state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor, time_list = compute_dimer_state_dilution_factor(file_path)

        dilution_factor = dilution_factor.tolist()
        Vt_list_states = [Vt] * state_num

        dilution_factor_list = dilution_factor_list + dilution_factor
        state_energy_list = state_energy_list + state_energy.tolist()
        Vt_for_state_list = Vt_for_state_list + Vt_list_states

    dilution_factor_list = np.array(dilution_factor_list)
    state_energy_list = np.array(state_energy_list)
    Vt_for_state_list = np.array(Vt_for_state_list)

    dilution_factor_criterion2 = 0.01
    length = len(dilution_factor_list)

    localization_index_list = [index for index in range(length) if
                               dilution_factor_list[index] > dilution_factor_criterion2]

    extended_state_index_list = [index for index in range(length) if
                                 dilution_factor_list[index] <= dilution_factor_criterion2]

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    ax.scatter(state_energy_list [localization_index_list] , Vt_for_state_list[localization_index_list] , marker = 'o' , color = 'red' , s=40)
    ax.scatter(state_energy_list[extended_state_index_list] , Vt_for_state_list[extended_state_index_list] , marker =  'v', color = 'black', s=40)
    # ax.scatter(T_list[mixing_index_list], T_prime_list[mixing_index_list], marker = 's' , color = 'purple', s = 40 )

    ax.set_xlabel("$E$")
    ax.set_ylabel(" $V_{t}$ ")

    if save_bool:
        fig_name = "dilution factor vs energy and Vt.svg"
        fig_name = os.path.join(folder_path, fig_name)
        fig.savefig(fig_name)

    plt.show()

