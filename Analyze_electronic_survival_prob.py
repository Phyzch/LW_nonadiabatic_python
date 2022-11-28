from util import *
from Analyze_survival_prob import Read_survival_prob_all_state, compute_dilution_factor
from anharmonic_transition_factor_T import estimate_anharmonic_transition_factor_T_state_lists
from nonadiabatic_transition_factor_T import analyze_nonadiabatic_transition_factor_lists

def Read_electronic_survival_prob_all_state(file_path):
    file_name = os.path.join(file_path,"electronic_survival_prob.txt")

    state_num = 0
    state_energy = []
    mode_number_list = []
    time_list = []
    survival_prob_list = []
    with open (file_name) as f:
        data = f.read()
        data = re.split('\n', data)
        data = [i for i in data if i != '']
        datalen = len(data)

        line = data[0]
        line = re.split(' ', line)
        state_num = int(line[0])

        line = data[1]
        line = re.split(' ' , line)
        state_energy = [float(i) for i in line if i!= '']

        line_index = 2
        for i in range(state_num):
            line = data[line_index]
            line = re.split(' ', line)
            mode_number = [int(i) for i in line if i!= '']

            mode_number_list.append(mode_number)
            line_index = line_index + 1

        step_num = int ((datalen - line_index)/2 )

        for i in range(step_num):
            line = data[line_index]
            line = re.split(' ', line)

            t = float(line[0])
            time_list.append(t)

            line_index = line_index  + 1

            line = data[line_index]
            line = re.split(' ', line)

            survival_prob = [float(i) for i in line if i!= '']

            survival_prob_list.append(survival_prob)

            line_index = line_index + 1

        state_energy = np.array(state_energy)
        mode_number_list = np.array(mode_number_list)
        survival_prob_list = np.array(survival_prob_list )
        time_list = np.array(time_list)

        nmode = len(mode_number_list[0])
        return state_num, nmode, state_energy, mode_number_list, survival_prob_list, time_list

def plot_electronic_survival_prob_with_vib_survival_prob():
    '''
    Compare survival probability
    :return:
    '''
    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl try/"

    state_num, nmode, state_energy, mode_number_list, survival_prob_list, time_list = Read_survival_prob_all_state(file_path)
    survival_prob_list = np.transpose(survival_prob_list)
    _,_,_,_, electronic_survival_prob_list , _ = Read_electronic_survival_prob_all_state(file_path)
    electronic_survival_prob_list = np.transpose(electronic_survival_prob_list)

    # take only vibrational mode number
    mode_number_list = mode_number_list[:,1:]

    color_list = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

    # plot survival prob and electronic survival prob. same state using same color
    fig = plt.figure(figsize = (20,10))
    spec = gridspec.GridSpec(nrows = 1, ncols = 2, figure = fig)
    ax1 = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])

    index_list = [0,2, 4, 6, 7]
    index_len = len(index_list)

    for i in range(index_len):
        index = index_list[i]

        ax1.plot(time_list, survival_prob_list[index], label = "n = " + str(mode_number_list[index]), color = color_list[i] , linewidth = 2)
        ax2.plot(time_list, electronic_survival_prob_list[index] , color = color_list[i] , linewidth = 2)

    ax1.set_yscale('log')
    ax1.set_ylabel('P(t)')
    ax1.set_xlabel('t')
    ax1.legend(loc = 'best')

    ax2.set_ylabel('electronic state $P(t)$')
    ax2.set_xlabel('t')

    plt.show()


def Analyze_electronic_dilution_factor_with_T_phase_diagram():
    '''
    See Leitner's review paper (2015) for estimation formula for T.
    :return:
    '''
    # anharmonic coupling
    V0 = 300
    anharmonic_scaling_factor = 0.3

    save_bool = True
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/" \
                  "Bchl 5mode/batch_simulation_phase_diagram/batch_simulation_all_with_electronic_survival_prob"

    file_path1 = "Vt=0"
    file_path2 = "Vt=100"
    file_path3 = "Vt=300"
    file_path4 = "Vt=500"

    file_path_list = [file_path1, file_path2, file_path3, file_path4]
    file_path_list = [os.path.join(folder_path, file_path) for file_path in file_path_list]

    Vt_list = [0,100,300,500]
    alpha_list = np.array([ 0.169, 0.163, 0.127, 0.101, 0.101])

    T_list = []
    T_prime_list = []
    electronic_dilution_factor_list = []

    path_num = len(file_path_list)

    for i in range(path_num):
        file_path = file_path_list[i]
        Vt = Vt_list[i]

        state_num, nmode, state_energy, mode_number_list, electronic_survival_prob, time_list = Read_electronic_survival_prob_all_state(file_path)
        electronic_survival_prob = np.transpose(electronic_survival_prob)
        # take average of survival prob at later time
        electronic_dilution_factor = np.mean(electronic_survival_prob[:, int(len(time_list)/2): ], 1)

        mode_number_list = [x[1:] for x in mode_number_list]

        # transition factor T in the same PES
        T = estimate_anharmonic_transition_factor_T_state_lists(V0, anharmonic_scaling_factor , mode_number_list)

        # transition factor in different PES
        T_prime = analyze_nonadiabatic_transition_factor_lists(mode_number_list, Vt, alpha_list)

        # concatenate lists
        T_list = T_list + T
        T_prime_list = T_prime_list + T_prime
        electronic_dilution_factor_list = electronic_dilution_factor_list + electronic_dilution_factor.tolist()

    electronic_dilution_factor_list = np.array(electronic_dilution_factor_list)
    T_list = np.array(T_list)
    T_prime_list = np.array(T_prime_list)

    electronic_dilution_factor_criterion1 = 0.8
    electronic_dilution_factor_criterion2 = 0.6
    length = len(electronic_dilution_factor_list)

    localization_index_list = [index for index in range(length) if
                               electronic_dilution_factor_list[index] > electronic_dilution_factor_criterion1]
    mixing_index_list = [index for index in range(length) if
                         electronic_dilution_factor_list[index] <= electronic_dilution_factor_criterion1
                         and electronic_dilution_factor_list[index] > electronic_dilution_factor_criterion2 ]

    extended_state_index_list = [index for index in range(length) if
                                 electronic_dilution_factor_list[index] <= electronic_dilution_factor_criterion2]

    # plot electronic dilution factor against T and T_prime (scatter plot)
    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    # gray color map
    # sc = ax.scatter( T_list, T_prime_list, marker = 'o' , s = 40, c = electronic_dilution_factor_list, cmap = 'gray_r' )
    # cbar = fig.colorbar(sc, ax = ax)
    # cbar.set_label("$\sigma_{e}$ ")

    # distinguish result with different color
    ax.scatter(T_list[localization_index_list] , T_prime_list[localization_index_list] , marker = 'o' , color = 'red' , s=40)
    ax.scatter(T_list[extended_state_index_list] , T_prime_list[extended_state_index_list] , marker =  'v', color = 'blue', s=40)
    ax.scatter(T_list[mixing_index_list], T_prime_list[mixing_index_list], marker = 's' , color = 'purple', s = 40 )


    T_boundary = np.linspace(0,1, 100)
    T_prime_boundary = 1 - T_boundary
    ax.plot(T_boundary, T_prime_boundary, color = 'black')
    ax.set_ylim([-0.05,0.5])
    ax.set_xlabel('T')
    ax.set_ylabel("T'")

    if save_bool:
        fig_name = "electronic dilution factor vs T.svg"
        fig_name = os.path.join(folder_path, fig_name)
        fig.savefig(fig_name)

    plt.show()