from Analyze_dimer_code.Analyze_dimer_batch_siimulation_survival_prob_auxiliary_func import *
from estimate_transition_factor_T_full_mode import estimate_transition_factor_with_freq_and_energy_dimer
from estimate_transition_factor_T_prime_full_mode import estimate_T_prime_prefactor_subroutine_dimer

def analyze_T_T_prime_vs_energy_main():
    '''

    :return:
    '''
    # Bigwood formula T
    # plot_anharmonic_T_estimate_for_data_vs_energy_V0_3050()

    # Bigwood formula T + T'
    plot_T_T_prime_estimate_for_data_vs_energy_V0_3050()
def plot_anharmonic_T_estimate_for_data_vs_energy_V0_3050():
    '''

    :return:
    '''
    save_bool = False
    fig_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/" \
               "5_mode/batch_simulation_Bigwood_scaling/batch_simulation_output/"
    V0 = 3050
    frequency_list = np.array([890, 727, 345, 1117, 1158])
    dof = len(frequency_list)

    scaling_factor_list = np.sqrt(frequency_list) / 270
    scaling_factor = np.prod(np.power(scaling_factor_list, 1 / dof))

    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/" \
                "5_mode/batch_simulation_Bigwood_scaling/batch_simulation_output/Vt=0/"



    # T_dimer = T_monomer1 + T_monomer2
    T_states = estimate_transition_factor_T_dimer_lists(V0, scaling_factor, file_path)

    # assume equal energy in each monomer.
    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor, time_list = compute_dimer_state_dilution_factor(
        file_path)
    max_energy = np.max(state_energy)
    min_energy = np.min(state_energy)
    energy_data_num = 40
    energy_list = np.linspace(min_energy, max_energy, energy_data_num)
    T_list_equal_energy_assumption = []
    for i in range(energy_data_num):
        energy = energy_list[i]
        T = estimate_transition_factor_with_freq_and_energy_dimer(energy, frequency_list, scaling_factor)
        T_list_equal_energy_assumption.append(T)


    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    ax.scatter(state_energy, T_states, marker='o', color='red', s=40)
    ax.plot(energy_list, T_list_equal_energy_assumption, color = 'black' , linewidth = 2 )

    ax.set_xlabel('$E$')
    ax.set_ylabel('$T$')

    if save_bool:
        fig_name = "anharmonic T vs energy E.svg"
        fig_name = os.path.join(fig_path, fig_name)
        fig.savefig(fig_name)

    plt.show()

def plot_T_T_prime_estimate_for_data_vs_energy_V0_3050():
    '''

    :return:
    '''
    save_bool = False
    fig_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/" \
               "5_mode/batch_simulation_Bigwood_scaling/batch_simulation_output/"

    frequency_list_monomer1 = np.array([890, 727, 345, 1117, 1158])
    frequency_list_monomer2 = np.array([890, 727, 345, 1117, 1158])

    monomer1_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])
    monomer2_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])

    dimer_alpha = [monomer1_alpha , monomer2_alpha]
    dimer_frequency = [ frequency_list_monomer1 , frequency_list_monomer2 ]

    V0 = 3050
    frequency_list = np.array([890, 727, 345, 1117, 1158])
    dof = len(frequency_list)

    scaling_factor_list = np.sqrt(frequency_list) / 270
    scaling_factor = np.prod(np.power(scaling_factor_list, 1 / dof))

    Vt = 363

    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/" \
                "5_mode/batch_simulation_Bigwood_scaling/batch_simulation_output/Vt=363/"

    # T_dimer = T_monomer1 + T_monomer2
    T_states = estimate_transition_factor_T_dimer_lists(V0, scaling_factor, file_path)
    T_prime_states = estimate_nonadiabatic_transition_factor_T_prime_dimer_list(dimer_alpha, dimer_frequency,
                                                                                file_path, Vt)

    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor, time_list = compute_dimer_state_dilution_factor(
        file_path)

    T_T_prime_sum_list = T_states + T_prime_states


    # estimate using equal energy per monomer assumption
    alpha_list = np.array([0.169, 0.163, 0.127, 0.101, 0.101])

    max_energy = np.max(state_energy)
    min_energy = np.min(state_energy)
    energy_data_num = 40
    energy_list = np.linspace(min_energy, max_energy, energy_data_num)
    T_list_equal_energy_assumption = []
    T_prime_list_equal_energy_assumption = []

    for i in range(energy_data_num):
        energy = energy_list[i]
        T = estimate_transition_factor_with_freq_and_energy_dimer(energy, frequency_list, scaling_factor)
        T_list_equal_energy_assumption.append(T)

        T_prime_prefactor = estimate_T_prime_prefactor_subroutine_dimer(alpha_list, frequency_list , energy)
        T_prime_estimate= T_prime_prefactor* Vt
        T_prime_list_equal_energy_assumption.append(T_prime_estimate)

    T_T_prime_sum_estimate = np.array(T_prime_list_equal_energy_assumption) + np.array(T_list_equal_energy_assumption)


    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    # T
    ax.scatter(state_energy, T_T_prime_sum_list , marker='o', color='red', s=40)
    ax.plot(energy_list, T_T_prime_sum_estimate, color = 'black')

    # T+ T'
    ax.scatter(state_energy, T_states, marker='o', color='blue', s=40)
    ax.plot(energy_list, T_list_equal_energy_assumption, color = 'black' , linewidth = 2 )

    ax.set_xlabel('$E$')
    ax.set_ylabel("transition factor")

    if save_bool:
        fig_name = "T+T' vs energy E.svg"
        fig_name = os.path.join(fig_path, fig_name)
        fig.savefig(fig_name)

    # figure out energy
    T_transition_energy_index_rank = np.argsort(np.abs(np.array(T_list_equal_energy_assumption)-1) )
    T_transition_energy = np.mean( energy_list[T_transition_energy_index_rank[:2]])

    T_T_prime_transition_energy_index_rank = np.argsort(np.abs(np.array(T_T_prime_sum_estimate)-1) )
    T_T_prime_transition_energy = np.mean( energy_list[ T_T_prime_transition_energy_index_rank[:2] ] )

    print("energy for T transition : " + str(T_transition_energy))
    print("energy for T' + T transition " + str(T_T_prime_transition_energy))

    plt.show()