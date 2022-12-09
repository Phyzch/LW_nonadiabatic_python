from util import *
from Analyze_dimer_code.Analyze_dimer_batch_siimulation_survival_prob_auxiliary_func import *

def analyze_nonstatistical_states():
    '''

    :return:
    '''

    find_out_localized_state_at_high_energy()

    # find_out_extended_state_at_low_energy()


def find_out_localized_state_at_high_energy():
    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/" \
                "batch_simulation_Bigwood_scaling/batch_simulation_output/Vt=0/"

    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor_all_states,\
        time_list = compute_dimer_state_dilution_factor(file_path)

    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, survival_prob_list,\
        time_list = Read_dimer_survival_prob_all_states(file_path)

    localized_state_high_energy = [index for index in range(state_num) if (dilution_factor_all_states[index] > 0.1 and
                                                                           state_energy[index] > 20000) ]

    # localized_state_high_energy = [index for index in range(state_num) if (dilution_factor_all_states[index] < 0.01 and
    #                                                                        state_energy[index] < 25000) ]

    print("index: " + str(localized_state_high_energy))

    for i in range(len(localized_state_high_energy)):
        monomer1 = monomer1_quantum_num[ localized_state_high_energy[i] ]
        monomer2 = monomer2_quantum_num [localized_state_high_energy[i]]
        es1, es2 = compute_edge_factor_ei(monomer1, monomer2)
        print( str(monomer1) + "   " + str(monomer2) )
        es_dimer = [round(es1,3),  round(es2,3)]
        print(es_dimer)

    fig = plt.figure(figsize = (10,10))
    spec = gridspec.GridSpec(nrows = 1, ncols = 1, figure = fig)
    ax = fig.add_subplot(spec[0,0])

    for i in range(len(localized_state_high_energy)):
        index = localized_state_high_energy[i]
        monomer1 = monomer1_quantum_num[index]
        monomer2 = monomer2_quantum_num[index]


        label = "|v1>=" + str(monomer1) + " |v2>=" + str(monomer2) + " E = " + str(state_energy[index])
        ax.plot(time_list, survival_prob_list[index] , label = label , linewidth = 2)


    # index = 13
    # label = "|v1>=" + str(monomer1_quantum_num[index]) + " |v2>=" + str(monomer2_quantum_num[index])
    # ax.plot(time_list, survival_prob_list[index], label = label , linewidth = 2)

    ax.legend(loc = 'best')
    ax.set_yscale('log')

    plt.show()

def find_out_extended_state_at_low_energy():
    '''

    :return:
    '''
    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/" \
                "batch_simulation_Bigwood_scaling/batch_simulation_output/Vt=363/"

    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor_all_states, \
        time_list = compute_dimer_state_dilution_factor(file_path)

    state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, survival_prob_list, \
        time_list = Read_dimer_survival_prob_all_states(file_path)

    exteneded_state_at_low_energy = [index for index in range(state_num) if (dilution_factor_all_states[index] < 0.02 and
                                                                           state_energy[index] < 10000)]

    for i in range(len(exteneded_state_at_low_energy)):
        monomer1 = monomer1_quantum_num[exteneded_state_at_low_energy[i]]
        monomer2 = monomer2_quantum_num[exteneded_state_at_low_energy[i]]
        print(str(monomer1) + "   " + str(monomer2))

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    # for i in range(len(exteneded_state_at_low_energy)):
    #     index = exteneded_state_at_low_energy[i]
    #     monomer1 = monomer1_quantum_num[index]
    #     monomer2 = monomer2_quantum_num[index]
    #
    #     label = "|v1>=" + str(monomer1) + " |v2>=" + str(monomer2)
    #     ax.plot(time_list, survival_prob_list[index], label=label, linewidth=2)


    index = 70
    label = "|v1>=" + str(monomer1_quantum_num[index]) + " |v2>=" + str(monomer2_quantum_num[index])
    ax.plot(time_list, survival_prob_list[index], label = label , linewidth = 2)

    ax.legend(loc='best')
    ax.set_yscale('log')

    plt.show()


def compute_edge_factor_ei(monomer1, monomer2):
    '''
     J. Chem. Phys. 130, 134310 (2009);  eq.6
    :param monomer1:
    :param monomer2:
    :return:
    '''
    N = len(monomer1)
    average_qn1 = np.mean(monomer1)
    es1 = np.sqrt( np.sum( np.power( monomer1 / average_qn1 - 1 , 2)) / (N * (N-1)) )
    average_qn2 = np.mean(monomer2)
    es2 = np.sqrt( np.sum( np.power( monomer2 / average_qn2 - 1 , 2)) / (N * (N-1)) )

    return es1, es2