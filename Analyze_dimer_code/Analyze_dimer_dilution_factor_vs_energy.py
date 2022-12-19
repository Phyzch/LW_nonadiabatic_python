from util import *

from dilution_factor_with_respect_to_barrier_height import compute_barrier_height
from Analyze_dimer_code.Analyze_dimer_batch_siimulation_survival_prob_auxiliary_func import *


def analyze_dilution_factor_vs_energy_main():
    '''

    :return:
    '''
    # V0 = 300, a = 0.3
    # analyze_dilution_factor_vs_energy_V0_300()

    # V0 = 3050, a = geometric mean , Bigwood formula (1998 Bigwood Gruebele Leitner Wolynes).
    # analyze_dilution_factor_vs_energy_V0_3050()

    # V0=3050, a = geometric mean. Bigwood formula (1998 Bigwood Gruebele Leitner Wolynes). energy in one monomer.
    analyze_dilution_factor_vs_energy_V0_3050_vib_energy_in_one_monomer()

def analyze_dilution_factor_vs_energy_V0_3050():
    '''

    :return:
    '''
    save_bool = False
    parent_folder = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/" \
                    "5_mode/batch_simulation_Bigwood_scaling/batch_simulation_output/"

    folder_path_no_coupling = "Vt=0"
    folder_path_no_coupling_list  = [folder_path_no_coupling]
    folder_path_no_coupling_list = [os.path.join(parent_folder, path) for path in folder_path_no_coupling_list]

    folder_path_with_Vt = "Vt=363"
    folder_path_with_Vt_list = [folder_path_with_Vt]
    folder_path_with_Vt_list = [os.path.join(parent_folder, path) for path in folder_path_with_Vt_list]

    Vt_list = [0, 363]

    # transition energy predicted by T=1
    transition_energy_T = 15800
    transition_energy_T_T_prime = 12400

    analyze_dilution_factor_vs_energy_subroutine(folder_path_no_coupling_list, folder_path_with_Vt_list, Vt_list,save_bool,
                                                 parent_folder, transition_energy_T,  transition_energy_T_T_prime)

def analyze_dilution_factor_vs_energy_V0_3050_vib_energy_in_one_monomer():
    '''

    :return:
    '''
    save_bool = False
    parent_folder = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/" \
                    "batch_simulation_Bigwood_scaling/batch_simulation_energy_in_one_monomer/output_file"

    folder_path_no_coupling = "Vt=0"
    folder_path_no_coupling_list  = [folder_path_no_coupling]
    folder_path_no_coupling_list = [os.path.join(parent_folder, path) for path in folder_path_no_coupling_list]

    folder_path_with_Vt = "Vt=363"
    folder_path_with_Vt_list = [folder_path_with_Vt]
    folder_path_with_Vt_list = [os.path.join(parent_folder, path) for path in folder_path_with_Vt_list]

    Vt_list = [0, 363]

    # transition energy predicted by T=1
    transition_energy_T = 12050
    transition_energy_T_T_prime = 9500

    analyze_dilution_factor_vs_energy_subroutine(folder_path_no_coupling_list, folder_path_with_Vt_list, Vt_list,save_bool,
                                                 parent_folder, transition_energy_T,  transition_energy_T_T_prime)

def analyze_dilution_factor_vs_energy_V0_300():
    '''

    :return:
    '''
    save_bool = False

    parent_folder = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/output_file/"

    folder_path_no_coupling1 = "Vt=0/"
    folder_path_no_coupling2 = "Vt=0_low_energy/"
    folder_path_no_coupling3 = "Vt=0_middle_energy/"
    folder_path_no_coupling4 = "Vt=0_high_energy/"
    folder_path_no_coupling_list  = [folder_path_no_coupling1 , folder_path_no_coupling2, folder_path_no_coupling3,
                                     folder_path_no_coupling4]

    folder_path_no_coupling_list = [os.path.join(parent_folder, path) for path in folder_path_no_coupling_list]

    folder_path_with_Vt1 = "Vt=363/"
    folder_path_with_Vt2 = "Vt=363_low_energy/"
    folder_path_with_Vt3 = "Vt=363_middle_energy/"
    folder_path_with_Vt4 = "Vt=363_high_energy/"

    folder_path_with_Vt_list = [folder_path_with_Vt1, folder_path_with_Vt2, folder_path_with_Vt3, folder_path_with_Vt4 ]
    folder_path_with_Vt_list = [os.path.join(parent_folder, path) for path in folder_path_with_Vt_list]

    Vt_list = [0, 363]

    analyze_dilution_factor_vs_energy_subroutine(folder_path_no_coupling_list, folder_path_with_Vt_list,
                                                 Vt_list, save_bool, parent_folder)


def analyze_dilution_factor_vs_energy_subroutine(folder_path_no_coupling_list, folder_path_with_Vt_list, Vt_list, save_bool, fig_folder, transition_energy_T,
                                                 transition_energy_T_T_prime):
    '''

    :return:
    '''
    path_num = len(folder_path_no_coupling_list)

    state_energy_list = []
    dilution_factor_no_coupling_list = []
    dilution_factor_with_Vt_list = []

    for j in range(2):
        if j == 0:
            path_list = folder_path_no_coupling_list
        else:
            path_list = folder_path_with_Vt_list

        for i in range(path_num):
            path = path_list[i]
            state_num, nmode, state_energy, monomer1_quantum_num, monomer2_quantum_num, dilution_factor, time_list = \
                compute_dimer_state_dilution_factor(path)

            state_energy = state_energy.tolist()
            dilution_factor = dilution_factor.tolist()
            if j == 0:
                dilution_factor_no_coupling_list = dilution_factor_no_coupling_list + dilution_factor
                state_energy_list = state_energy_list + state_energy
            else:
                dilution_factor_with_Vt_list = dilution_factor_with_Vt_list + dilution_factor

    state_energy_list = np.array(state_energy_list)
    dilution_factor_no_coupling_list = np.array(dilution_factor_no_coupling_list)
    dilution_factor_with_Vt_list = np.array(dilution_factor_with_Vt_list)


    marker_list = ["^", "s"]
    # plot dilution factor without coupling black , dilution factor with coupling red
    color_list = ['black', 'red']

    fig = plt.figure(figsize=(20, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=2, figure=fig)
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[0, 1])

    ax1.scatter(state_energy_list, dilution_factor_no_coupling_list, color = color_list[0], marker = marker_list[0] , label = "Vt = " + str(Vt_list[0]) )
    ax1.scatter(state_energy_list, dilution_factor_with_Vt_list, color = color_list[1] , marker = marker_list[1] , label = "Vt = " + str(Vt_list[1]))

    ax1.set_yscale('log')
    ax1.legend(loc='best', prop={'size': 15})

    ax1.set_xlabel('E')
    ax1.set_ylabel('$\sigma$')

    # indicate transition energy.
    ax1.axvline(x = transition_energy_T, linewidth = 3, color = 'black')
    ax1.axvline(x = transition_energy_T_T_prime , linewidth = 3 , color = 'brown')

    dilution_factor_ratio = dilution_factor_with_Vt_list / dilution_factor_no_coupling_list
    ax2.scatter(state_energy_list, dilution_factor_ratio, color = 'black' )

    ax2.set_xlabel('E')
    ax2.set_ylabel('$\sigma$($V_{t}=$' + str(Vt_list[1]) + ") / $\sigma$($V_{t}=0$)")
    ax2.legend(loc='best', prop={'size': 15})
    ax2.set_yscale('log')


    plt.show()

    if save_bool:
        fig_name = "dilution factor vs energy.svg"
        fig_path = os.path.join(fig_folder, fig_name)
        fig.savefig(fig_path)


