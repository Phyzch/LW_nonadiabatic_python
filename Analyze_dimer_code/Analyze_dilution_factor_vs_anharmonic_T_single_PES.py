from Analyze_dimer_code.Analyze_dimer_batch_siimulation_survival_prob_auxiliary_func import \
    compute_dimer_state_dilution_factor, estimate_transition_factor_T_dimer_lists
from util import *


def analyze_dilution_factor_with_anharmonic_T_dimer_main():
    '''

    :return:
    '''
    # V0 = 300 , a = 0.3
    # analyze_dilution_factor_with_anharmonic_T_V0_300_dimer()

    # V0 = 3050, a = geometric mean of 'a' for each mode
    # analyze_dilution_factor_with_anharmonic_T_V0_3050_dimer()

    # V0 = 3050, a = geometric mean of 'a' for each mode. energy in single monomer.
    analyze_dilution_factor_with_anharmonic_T_V0_3050_dimer_energy_in_one_monomer()


def analyze_dilution_factor_with_anharmonic_T_V0_3050_dimer():
    '''

    :return:
    '''
    V0 = 3050
    frequency_list =  np.array([890, 727, 345, 1117, 1158])
    dof = len(frequency_list)

    scaling_factor_list = np.sqrt(frequency_list) / 270
    scaling_factor = np.prod( np.power(scaling_factor_list , 1/dof) )

    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation_Bigwood_scaling/batch_simulation_output/"
    file_path = "Vt=0"
    file_path_list = [file_path]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    analyze_dilution_factor_with_anharmonic_T_dimer_subroutine(file_path_list, V0, scaling_factor, save_bool, folder_path)

def  analyze_dilution_factor_with_anharmonic_T_V0_3050_dimer_energy_in_one_monomer():
    '''

    :return:
    '''
    V0 = 3050
    frequency_list =  np.array([890, 727, 345, 1117, 1158])
    dof = len(frequency_list)

    scaling_factor_list = np.sqrt(frequency_list) / 270
    scaling_factor = np.prod( np.power(scaling_factor_list , 1/dof) )

    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/" \
                  "batch_simulation_Bigwood_scaling/batch_simulation_energy_in_one_monomer/output_file"

    file_path = "Vt=0_high_energy_dE=1500"
    file_path_list = [file_path]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    analyze_dilution_factor_with_anharmonic_T_dimer_subroutine(file_path_list, V0, scaling_factor, save_bool, folder_path)

def analyze_dilution_factor_with_anharmonic_T_V0_300_dimer():
    '''

    :return:
    '''
    V0 = 300
    scaling_factor = 0.3

    save_bool = True
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/output_file/"

    file_path1 = "Vt=0"
    file_path2 = "Vt=0_low_energy"
    file_path3 = "Vt=0_middle_energy"

    file_path_list = [file_path1, file_path2, file_path3]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    analyze_dilution_factor_with_anharmonic_T_dimer_subroutine(file_path_list, V0, scaling_factor, save_bool, folder_path)


def analyze_dilution_factor_with_anharmonic_T_dimer_subroutine(file_path_list, V0, scaling_factor, save_bool, fig_path ):
    T_list = []
    dilution_factor_list = []
    file_num = len(file_path_list)

    for i in range(file_num):
        file_path = file_path_list[i]
        state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor, time_list = compute_dimer_state_dilution_factor(
            file_path)

        dilution_factor = dilution_factor.tolist()

        # transition factor in same PES
        T = estimate_transition_factor_T_dimer_lists(V0, scaling_factor, file_path)
        T = T.tolist()

        T_list = T_list + T
        dilution_factor_list = dilution_factor_list + dilution_factor

    T_list = np.array(T_list)
    dilution_factor_list = np.array(dilution_factor_list)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    ax.scatter(T_list, dilution_factor_list, marker='o', s=40)

    ax.set_xlabel('T')
    ax.set_ylabel('$\sigma$')
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim([0.1 , np.max(T_list) + 0.5])
    ax.set_ylim([np.min(dilution_factor_list) / 1.5 , 1.2])

    plt.show()

    if save_bool:
        fig_name = "anharmonic T vs dilution factor.svg"
        fig_name = os.path.join(fig_path, fig_name)
        fig.savefig(fig_name)
