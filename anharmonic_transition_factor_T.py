from util import *
from Analyze_Nloc_with_dilution_factor import read_Nloc_anharmonic_nonadiabatic_coupling

def estimate_anharmonic_transition_factor_T_order_Q( Q, V0, scaling_factor , state_qn):
    '''
    formula we use: Leitner & Wolynes 1996 https://aip.scitation.org/doi/10.1063/1.472920.  eq.(4.6)
    See also section 2.3.2 in Leitner 2015: Quantum ergodicity and energy flow in molecules.
    Q : order of coupling. Usually Q=3,4 is largest contribution
    :return:
    '''
    frequency_list = np.array([ 890, 727, 345, 1117, 1158 ])
    omega_rms = np.sqrt( np.mean( np.power(frequency_list,2) ) )
    dof = len(frequency_list)
    # connectivity
    K = np.power( 2 * dof, Q) / np.math.factorial(Q)

    # local density of states
    Dq = 1/(np.pi * np.sqrt(Q) * omega_rms)

    # coupling strength
    energy = np.sum( frequency_list * np.array(state_qn) )
    M = energy / (dof * omega_rms)

    ladder_operator_contribution = np.power(M , Q/2)
    Vq = V0 * np.power(scaling_factor , Q) * ladder_operator_contribution

    Tq = np.sqrt(2 * np.pi / 3) * K * Dq * Vq

    return Tq

def estimate_anharmonic_transition_factor_T( V0, scaling_factor , state_qn):
    '''

    :param V0:
    :param scaling_factor:
    :param state_qn:
    :return:
    '''
    Q = 3
    T3 = estimate_anharmonic_transition_factor_T_order_Q(Q, V0, scaling_factor, state_qn)

    # Q = 4
    # T4 = estimate_anharmonic_transition_factor_T_order_Q(Q, V0, scaling_factor, state_qn)

    # for Q=4, estimation for K is problematic, therefore, we only take Q=3
    T = T3

    return T

def estimate_anharmonic_transition_factor_T_method2(Vq):
    '''

    :param V:
    :param state_qn:
    :return:
    '''
    Q  = 3 # order or coupling
    frequency_list = np.array([ 890, 727, 345, 1117, 1158 ])
    omega_rms = np.sqrt( np.mean( np.power(frequency_list,2) ) )
    dof = len(frequency_list)

    # connectivity
    K = np.power(2 * dof, Q) / np.math.factorial(Q)

    # local density of states
    Dq = 1/(np.pi * np.sqrt(Q) * omega_rms)

    Tq = np.sqrt(2 * np.pi / 3) * K * Dq * Vq

    return Tq

def estimate_anharmonic_transition_factor_T_single_state():
    '''

    :return:
    '''
    state = [0,2,2,2,0]

    V0 = 300
    scaling_factor = 0.3

    Q = 3
    T3 = estimate_anharmonic_transition_factor_T_order_Q(Q, V0, scaling_factor, state)

    # Q = 4
    # T4 = estimate_anharmonic_transition_factor_T_order_Q(Q, V0, scaling_factor, state)  # Q=4, estimation for K is problematic.

    T = T3

    # print ("T3 : " + str(T3) + " T4: "+ str(T4) )
    print(" T:" + str(T))

def compare_anharmonic_Nloc_and_transition_factor_T():
    '''

    :return:
    '''
    save_bool = False
    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_phase_diagram/states Nloc vs T data/"
    # file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl try/"

    V0 = 300
    scaling_factor = 0.3


    # read Nloc
    state_num, mode_num, mode_index_list, same_PES_Nloc_list, diff_PES_Nloc = read_Nloc_anharmonic_nonadiabatic_coupling(
        file_path)

    # transition factor T
    T_list = []
    for j in range(state_num):
        qn = mode_index_list[j]
        qn = qn[1:]
        T_ele = estimate_anharmonic_transition_factor_T(V0, scaling_factor, qn)
        T_list.append(T_ele)

    T_list_log = np.log(T_list)
    Nloc_log = np.log(same_PES_Nloc_list)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])


    ax.scatter( T_list_log , Nloc_log , edgecolors = 'black', facecolors = 'none', s=50 )
    ax.set_xlabel('log(T)')
    ax.set_ylabel('log($N_{loc}$)')

    # ax.scatter(T_list, same_PES_Nloc_list,  edgecolors = 'black', facecolors = 'none', s=50)
    # ax.set_xlabel('T')
    # ax.set_ylabel('$N_{loc}$')

    if save_bool:
        fig_name = "anharmonic Nloc vs T.svg"
        fig_name = os.path.join(file_path, fig_name)
        fig.savefig(fig_name)

    plt.show()

# estimate_anharmonic_transition_factor_T_single_state()