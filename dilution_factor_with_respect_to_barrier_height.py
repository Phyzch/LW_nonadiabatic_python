from util import *
from Analyze_survival_prob import compute_dilution_factor

def compute_barrier_height(ground_state_energy_shift , EV_coupling_alpha):
    '''
    # ground state energy for electronic state |1> compared to electronic state 0.
    # electron vibration coupling factor: EV_coupling_alpha =
    :return:
    '''
    frequency_list = np.array( [890, 727, 345, 1117, 1158] )

    # ci : coupling strength
    EV_coupling = EV_coupling_alpha * np.power(frequency_list , 3/2) * np.sqrt(2)

    # find ground state energy on the crossing surface.
    # fraction: c_{i}^{2}/(m omega_{i}^{2})
    strength_factor = np.power(EV_coupling , 2) / np.power(frequency_list,2)
    strength_fraction = strength_factor / np.sum(strength_factor)

    energy_barrier_height = np.sum( 0.5 * np.power(frequency_list, 2) * np.power( EV_coupling/(2 * np.power(frequency_list,2)) + 1/EV_coupling * strength_fraction * ground_state_energy_shift , 2 ) )

    return energy_barrier_height


def plot_dilution_factor_vs_energy():
    '''

    :return:
    '''

    save_bool = False

    ground_state_energy_shift = 600
    EV_coupling_alpha = np.array([0.169, 0.163, 0.127, 0.101, 0.101])
    energy_barrier_height = compute_barrier_height(ground_state_energy_shift, EV_coupling_alpha)

    parent_folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_phase_diagram/"

    folder_path_no_coupling1 = "batch_simulation_dE=600/Vt=0,V0=300,a=0.3/"
    folder_path_no_coupling2 = "batch_simulation_dE=600_high_energy/Vt=0/"
    folder_path_no_coupling3 = "batch_simulation_dE=600_low_energy/Vt=0/"
    folder_path_no_coupling_list = [folder_path_no_coupling1 , folder_path_no_coupling2, folder_path_no_coupling3]
    folder_path_no_coupling_list = [ os.path.join(parent_folder_path, path) for path in folder_path_no_coupling_list ]

    folder_path_with_Vt1 = "batch_simulation_dE=600/Vt=300,V0=300,a=0.3/"
    folder_path_with_Vt2 = "batch_simulation_dE=600_high_energy/Vt=300/"
    folder_path_with_Vt3 = "batch_simulation_dE=600_low_energy/Vt=300/"

    folder_path_with_Vt_list = [folder_path_with_Vt1 , folder_path_with_Vt2, folder_path_with_Vt3]
    folder_path_with_Vt_list = [os.path.join(parent_folder_path, path) for path in folder_path_with_Vt_list]

    Vt_list = [0,300]

    path_num = len(folder_path_no_coupling_list)

    state_energy_list = []
    electronic_state_list = []
    dilution_factor_no_coupling_list = []
    dilution_factor_with_Vt_list = []

    for j in range(2):
        if j == 0:
            path_list = folder_path_no_coupling_list
        else:
            path_list = folder_path_with_Vt_list

        for i in range(path_num):
            path = path_list[i]

            vib_state_energy, mode_number, dilution_factor = compute_dilution_factor(path)

            electronic_state = np.array( [ mode_number_ele[0] for mode_number_ele in mode_number ] )
            state_energy = vib_state_energy + electronic_state * ground_state_energy_shift

            if j == 0:
                electronic_state_list = electronic_state_list + electronic_state.tolist()
                state_energy_list = state_energy_list + state_energy.tolist()

            if j == 0:
                dilution_factor_no_coupling_list = dilution_factor_no_coupling_list + dilution_factor.tolist()
            else:
                dilution_factor_with_Vt_list = dilution_factor_with_Vt_list + dilution_factor.tolist()

    state_energy_list = np.array(state_energy_list)
    electronic_state_list = np.array(electronic_state_list)
    dilution_factor_no_coupling_list = np.array(dilution_factor_no_coupling_list)
    dilution_factor_with_Vt_list = np.array(dilution_factor_with_Vt_list)

    # plot state on electronic state 1 with triangle symbol , state on electronic state 2 with square symbol
    marker_list = ["^" , "s"]
    # plot dilution factor without coupling black , dilution factor with coupling red
    color_list = ['black' , 'red']

    electronic0_state_list = [index for index in range(len(electronic_state_list)) if electronic_state_list[index] == 0]
    electronic1_state_list = [index for index in range(len(electronic_state_list)) if electronic_state_list[index] == 1]

    fig = plt.figure(figsize=(20, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=2, figure=fig)
    ax1 = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])

    ax1.scatter( state_energy_list[electronic0_state_list], dilution_factor_no_coupling_list[electronic0_state_list] ,
                color= color_list[0] , marker = marker_list[0] , label = "Vt = " + str(Vt_list[0]) + " electronic state |0>" )
    ax1.scatter( state_energy_list[electronic1_state_list] , dilution_factor_no_coupling_list[electronic1_state_list] ,
                color= color_list[0], marker = marker_list[1], label = "Vt = " + str(Vt_list[0]) + " electronic state |1>" )

    ax1.scatter(state_energy_list[electronic0_state_list], dilution_factor_with_Vt_list[electronic0_state_list],
               color=color_list[1], marker=marker_list[0] , label = "Vt = " + str(Vt_list[1]) + " electronic state |0>")
    ax1.scatter(state_energy_list[electronic1_state_list], dilution_factor_with_Vt_list[electronic1_state_list],
               color=color_list[1], marker=marker_list[1] , label = "Vt = " + str(Vt_list[1]) + " electronic state |1>")

    ax1.axvline( x = energy_barrier_height , linewidth = 3)
    ax1.set_yscale('log')
    ax1.legend(loc = 'best' , prop={'size': 10})

    ax1.set_xlabel('E')
    ax1.set_ylabel('$\sigma$')

    # plot decrease of dilution factor by ratio in ax2
    dilution_factor_ratio = dilution_factor_with_Vt_list / dilution_factor_no_coupling_list

    ax2.scatter( state_energy_list[electronic0_state_list], dilution_factor_ratio[electronic0_state_list],
                 color = 'black' , marker = marker_list[0] , label = 'electronic state 0' )

    ax2.scatter( state_energy_list[electronic1_state_list], dilution_factor_ratio[electronic1_state_list],
                 color = 'red' , marker = marker_list[1], label = 'electronic state 1' )

    ax2.set_xlabel('E')
    ax2.set_ylabel('$\sigma$($V_{t}=$' + str(Vt_list[1]) + ") / $\sigma$($V_{t}=0$)" )
    ax2.axvline( x = energy_barrier_height, linewidth = 3)
    ax2.legend(loc = 'best', prop={'size': 10})

    plt.show()

    if save_bool:
        fig_name = "dilution factor vs energy.svg"
        fig_path = os.path.join(parent_folder_path , fig_name)
        fig.savefig(fig_path)

def plot_dilution_factor_vs_energy_strong_EV_coupling():
    '''

    :return:
    '''

    save_bool = True

    ground_state_energy_shift = 600
    EV_coupling_alpha = np.array([0.4,0.4,0.4,0.4,0.4])
    energy_barrier_height = compute_barrier_height(ground_state_energy_shift , EV_coupling_alpha)

    parent_folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/strong EV coupling model/large_qn_space/"

    folder_path_no_coupling1 = "Vt=0/"
    folder_path_no_coupling_list = [folder_path_no_coupling1]
    folder_path_no_coupling_list = [ os.path.join(parent_folder_path, path) for path in folder_path_no_coupling_list ]

    folder_path_with_Vt1 = "Vt=200/"

    folder_path_with_Vt_list = [folder_path_with_Vt1]
    folder_path_with_Vt_list = [os.path.join(parent_folder_path, path) for path in folder_path_with_Vt_list]

    Vt_list = [0,200]

    path_num = len(folder_path_no_coupling_list)

    state_energy_list = []
    electronic_state_list = []
    dilution_factor_no_coupling_list = []
    dilution_factor_with_Vt_list = []

    for j in range(2):
        if j == 0:
            path_list = folder_path_no_coupling_list
        else:
            path_list = folder_path_with_Vt_list

        for i in range(path_num):
            path = path_list[i]

            vib_state_energy, mode_number, dilution_factor = compute_dilution_factor(path)

            electronic_state = np.array( [ mode_number_ele[0] for mode_number_ele in mode_number ] )
            state_energy = vib_state_energy + electronic_state * ground_state_energy_shift

            if j == 0:
                electronic_state_list = electronic_state_list + electronic_state.tolist()
                state_energy_list = state_energy_list + state_energy.tolist()

            if j == 0:
                dilution_factor_no_coupling_list = dilution_factor_no_coupling_list + dilution_factor.tolist()
            else:
                dilution_factor_with_Vt_list = dilution_factor_with_Vt_list + dilution_factor.tolist()

    state_energy_list = np.array(state_energy_list)
    electronic_state_list = np.array(electronic_state_list)
    dilution_factor_no_coupling_list = np.array(dilution_factor_no_coupling_list)
    dilution_factor_with_Vt_list = np.array(dilution_factor_with_Vt_list)

    # plot state on electronic state 1 with triangle symbol , state on electronic state 2 with square symbol
    marker_list = ["^" , "s"]
    # plot dilution factor without coupling black , dilution factor with coupling red
    color_list = ['black' , 'red']

    electronic0_state_list = [index for index in range(len(electronic_state_list)) if electronic_state_list[index] == 0]
    electronic1_state_list = [index for index in range(len(electronic_state_list)) if electronic_state_list[index] == 1]

    fig = plt.figure(figsize=(20, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=2, figure=fig)
    ax1 = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])

    ax1.scatter( state_energy_list[electronic0_state_list], dilution_factor_no_coupling_list[electronic0_state_list] ,
                color= color_list[0] , marker = marker_list[0] , label = "Vt = " + str(Vt_list[0]) + " electronic state |0>" )
    ax1.scatter( state_energy_list[electronic1_state_list] , dilution_factor_no_coupling_list[electronic1_state_list] ,
                color= color_list[0], marker = marker_list[1], label = "Vt = " + str(Vt_list[0]) + " electronic state |1>" )

    ax1.scatter(state_energy_list[electronic0_state_list], dilution_factor_with_Vt_list[electronic0_state_list],
               color=color_list[1], marker=marker_list[0] , label = "Vt = " + str(Vt_list[1]) + " electronic state |0>")
    ax1.scatter(state_energy_list[electronic1_state_list], dilution_factor_with_Vt_list[electronic1_state_list],
               color=color_list[1], marker=marker_list[1] , label = "Vt = " + str(Vt_list[1]) + " electronic state |1>")

    ax1.axvline( x = energy_barrier_height , linewidth = 3)
    ax1.set_yscale('log')
    ax1.legend(loc = 'best' , prop={'size': 10})

    ax1.set_xlabel('E')
    ax1.set_ylabel('$\sigma$')

    # plot decrease of dilution factor by ratio in ax2
    dilution_factor_ratio = dilution_factor_with_Vt_list / dilution_factor_no_coupling_list

    ax2.scatter( state_energy_list[electronic0_state_list], dilution_factor_ratio[electronic0_state_list],
                 color = 'black' , marker = marker_list[0] , label = 'electronic state 0' )

    ax2.scatter( state_energy_list[electronic1_state_list], dilution_factor_ratio[electronic1_state_list],
                 color = 'red' , marker = marker_list[1], label = 'electronic state 1' )

    ax2.set_xlabel('E')
    ax2.set_ylabel('$\sigma$($V_{t}=$' + str(Vt_list[1]) + ") / $\sigma$($V_{t}=0$)" )
    ax2.axvline( x = energy_barrier_height, linewidth = 3)
    ax2.legend(loc = 'best', prop={'size': 10})

    plt.show()

    if save_bool:
        fig_name = "dilution factor vs energy.svg"
        fig_path = os.path.join(parent_folder_path , fig_name)
        fig.savefig(fig_path)