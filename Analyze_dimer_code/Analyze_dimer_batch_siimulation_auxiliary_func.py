from util import *
from Analyze_survival_prob import Read_survival_prob_all_state
from anharmonic_transition_factor_T import estimate_anharmonic_transition_factor_T
from Overlap_of_displaced_state import effective_num_coupling_submodule

def Read_dimer_survival_prob_all_states(file_path):
    state_num, nmode, state_energy, mode_number_list, survival_prob_list, time_list = Read_survival_prob_all_state(file_path)

    # [state_num, time_len]
    survival_prob_list = np.transpose(survival_prob_list)

    nmode = int( nmode / 2 )
    monomer1_quantum_num = mode_number_list[:, : nmode]
    monomer2_quantum_num = mode_number_list[:, nmode:]

    return state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, survival_prob_list, time_list


def compute_dimer_state_dilution_factor(file_path):
    '''
    Average over survival probability at late time ( t > final_time / 2) to get survival probability
    :param file_path:
    :return:
    '''
    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, survival_prob_list, time_list = Read_dimer_survival_prob_all_states(file_path)

    time_len = len(time_list)
    dilution_factor_all_states = np.mean( survival_prob_list[:, int(time_len/2) : ] , 1 )

    return state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor_all_states, time_list


def estimate_transition_factor_T_dimer_lists(V0, scaling_factor, file_path):
    '''

    :return:
    '''
    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, survival_prob_list, time_list = Read_dimer_survival_prob_all_states(file_path)
    T_list = estimate_transition_factor_T_dimer(V0, scaling_factor, state_num, monomer1_quantum_num, monomer2_quantum_num)

    return T_list

def estimate_transition_factor_T_dimer(V0, scaling_factor, state_num, monomer1_quantum_num, monomer2_quantum_num):
    '''

    :return:
    '''
    T_list = []

    for i in range(state_num):
        qn_monomer1 = monomer1_quantum_num[i]
        qn_monomer2 = monomer2_quantum_num[i]

        T1 = estimate_anharmonic_transition_factor_T(V0,scaling_factor, qn_monomer1)
        T2 = estimate_anharmonic_transition_factor_T(V0,scaling_factor, qn_monomer2)

        T = T1 + T2

        T_list.append(T)

    T_list = np.array(T_list)

    return T_list

def analyze_nonadiabatic_transition_factor_dimer(monomer1_qn, monomer2_qn, monomer1_alpha, monomer2_alpha, frequency_monomer1, frequency_monomer2, Vt):
    '''
    Vt: nonadiabatic coupling strength
    qn_list : quantum number
    :return:
    '''
    mode_num1 = len(monomer1_qn)
    effective_num_list_monomer1 = np.zeros([mode_num1])
    for i in range(mode_num1):
        alpha = monomer1_alpha[i]
        qn = monomer1_qn[i]
        effective_num = effective_num_coupling_submodule(alpha, qn)
        effective_num_list_monomer1[i] = effective_num

    effective_num_prod_monomer1 = np.prod(effective_num_list_monomer1)

    mode_num2 = len(monomer2_qn)
    effective_num_list_monomer2 = np.zeros([mode_num2])
    for i in range(mode_num2):
        alpha = monomer2_alpha[i]
        qn = monomer2_qn[i]
        effective_num = effective_num_coupling_submodule(alpha, qn)
        effective_num_list_monomer2[i] = effective_num

    effective_num_prod_monomer2 = np.prod(effective_num_list_monomer2)

    # compute nonadiabatic transition factor T
    freq_monomer1_rms = np.sqrt( np.mean( np.power(frequency_monomer1 , 2) ))
    freq_monomer2_rms = np.sqrt( np.mean( np.power(frequency_monomer2 , 2) ))
    D = 1 / (np.pi *  freq_monomer1_rms)
    K = effective_num_prod_monomer1 * effective_num_prod_monomer2

    # Franck condon factor is approximated as 1/sqrt(K)
    T_prime = np.sqrt(2 * np.pi / 3) * D * Vt * np.sqrt(K)

    return T_prime
def estimate_nonadiabatic_transition_factor_T_prime_dimer( monomer1_alpha, monomer2_alpha, frequency_list_monomer1, frequency_list_monomer2,  state_num, monomer1_quantum_num, monomer2_quantum_num, Vt):
    '''

    :return:
    '''

    T_prime_list = []
    for i in range(state_num):
        qn_monomer1 = monomer1_quantum_num[i]
        qn_monomer2 = monomer2_quantum_num[i]
        T_prime = analyze_nonadiabatic_transition_factor_dimer(qn_monomer1, qn_monomer2, monomer1_alpha, monomer2_alpha, frequency_list_monomer1, frequency_list_monomer2, Vt )
        T_prime_list.append(T_prime)

    T_prime_list = np.array(T_prime_list)

    return T_prime_list

def estimate_nonadiabatic_transition_factor_T_prime_dimer_list(dimer_alpha, dimer_frequency, file_path, Vt):
    '''

    :param dimer_alpha:
    :param dimer_frequency:
    :param file_path:
    :return:
    '''
    monomer1_alpha = dimer_alpha[0]
    monomer2_alpha = dimer_alpha[1]

    monomer1_frequency = dimer_frequency[0]
    monomer2_frequency = dimer_frequency[1]

    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, survival_prob_list, time_list = Read_dimer_survival_prob_all_states(file_path)

    T_prime_list = estimate_nonadiabatic_transition_factor_T_prime_dimer(monomer1_alpha, monomer2_alpha,
                                                                         monomer1_frequency, monomer2_frequency,
                                                                         state_num, monomer1_quantum_num, monomer2_quantum_num, Vt )

    return T_prime_list