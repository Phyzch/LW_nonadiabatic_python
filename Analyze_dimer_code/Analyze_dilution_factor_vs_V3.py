'''
Analyze effect of dilution factor vs V3 and a
'''
from util import *
from Analyze_dimer_code.Analyze_dimer_batch_siimulation_survival_prob_auxiliary_func import *

def dilution_factor_vs_anharmonic_coupling_nonadiabatic_system_main():
    '''

    :return:
    '''

    # plot_dilution_factor_vs_V3()

    plot_dilution_factor_vs_scaling_factor_a()

def plot_dilution_factor_vs_V3():
    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/test_anharmonic_effect/dilution_factor_vs_V3/output/"

    file_path1 = "V0=0"
    file_path2 = "V0=100,a=0.25"
    file_path3 = "V0=150,a=0.25"
    file_path4 = "V0=200,a=0.25"
    file_path5 = "V0=250,a=0.25"
    file_path6 = "V0=300,a=0.25"
    file_path7 = "V0=400,a=0.25"

    path_list = [file_path1, file_path2, file_path3, file_path4, file_path5, file_path6, file_path7]
    path_list = [os.path.join(folder_path, path) for path in path_list]

    V0_list = [0, 100, 150, 200, 250, 300, 400]
    scaling_factor = 0.25

    V3_list = np.array(V0_list) * np.power(scaling_factor,3)

    path_num = len(path_list)

    dilution_factor_all_states = []
    for i in range(path_num):
        file_path = path_list[i]

        state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, \
        dilution_factor, time_list = compute_dimer_state_dilution_factor(file_path)

        dilution_factor_all_states.append(dilution_factor)

    dilution_factor_all_states = np.array(dilution_factor_all_states)
    dilution_factor_all_states = np.transpose(dilution_factor_all_states)

    state_num = len(monomer1_quantum_num)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    for state_index in range(state_num):
        dilution_factor = dilution_factor_all_states[state_index]
        monomer1_qn = monomer1_quantum_num[state_index]
        monomer2_qn = monomer2_quantum_num[state_index]
        label = "v1=" + str(monomer1_qn) + "  v2= " + str(monomer2_qn)
        ax.plot(V3_list, dilution_factor, marker = 'o' , label = label, linewidth = 2, markersize = 10)

    ax.set_yscale('log')
    ax.legend(loc = 'best')
    ax.set_ylim([pow(10,-3) , 1])

    ax.set_xlabel('$V_{3}$  $cm^{-1}$')
    ax.set_ylabel('$\sigma$')

    if save_bool:
        fig_name = "dilution factor vs V3.svg"
        fig_path = os.path.join(folder_path, fig_name)
        fig.savefig(fig_path)

def plot_dilution_factor_vs_scaling_factor_a():
    '''

    :return:
    '''
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/" \
                "5_mode/batch_simulation/test_anharmonic_effect/dilution_factor_vs_a/output/"

    save_bool = False

    file_path1 = "a=0"
    file_path2 = "a=0.05"
    file_path3 = "a=0.1"
    file_path4 = "a=0.15"
    file_path5 = "a=0.2"
    file_path6 = "a=0.25"
    file_path7 = "a=0.3"

    file_path_list = [file_path1, file_path2, file_path3, file_path4, file_path5, file_path6, file_path7]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    V0 = 300
    scaling_factor_list = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]

    file_num = len(file_path_list)
    dilution_factor_all_states = []
    for i in range(file_num):
        file_path = file_path_list[i]

        state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, \
            dilution_factor, time_list = compute_dimer_state_dilution_factor(file_path)

        dilution_factor_all_states.append(dilution_factor)

    dilution_factor_all_states = np.array(dilution_factor_all_states)
    dilution_factor_all_states = np.transpose(dilution_factor_all_states)

    state_num = len(monomer1_quantum_num)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    for state_index in range(state_num):
        dilution_factor = dilution_factor_all_states[state_index]
        monomer1_qn = monomer1_quantum_num[state_index]
        monomer2_qn = monomer2_quantum_num[state_index]
        label = "v1=" + str(monomer1_qn) + "  v2= " + str(monomer2_qn)
        ax.plot(scaling_factor_list, dilution_factor, marker = 'o' , label = label, linewidth = 2, markersize = 10)

    ax.set_yscale('log')
    ax.legend(loc = 'best')
    ax.set_ylim([pow(10,-3) , 1])

    ax.set_xlabel('$a$ ')
    ax.set_ylabel('$\sigma$')

    if save_bool:
        fig_name = "dilution factor vs scaling factor.svg"
        fig_path = os.path.join(folder_path, fig_name)
        fig.savefig(fig_path)
