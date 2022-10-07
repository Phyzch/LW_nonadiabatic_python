from util import *

def Read_survival_prob_all_state(file_path):
    file_name = os.path.join(file_path,"survival_prob.txt")

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

def Analyze_survival_prob_all_state(file_path , save_bool):
    matplotlib.rcParams.update({'font.size': 15})
    state_num, nmode , state_energy, mode_number_list, survival_prob_list , time_list = Read_survival_prob_all_state(file_path)

    ground_state_index = [ i for i in range(state_num) if mode_number_list[i][0] == 0 ]
    excited_state_index = [i for i in range(state_num) if mode_number_list[i][0] == 1 ]

    State_energy_for_ground_state = state_energy[ground_state_index]
    State_energy_for_excited_state = state_energy[excited_state_index]

    # take IPR at final time as final IPR
    time_len = len(time_list)
    start_time_to_plot_index = int(time_len * 2 / 3)

    final_survival_prob = np.mean(survival_prob_list[ start_time_to_plot_index: ] , 0)
    final_survival_prob_ground_state = final_survival_prob[ground_state_index]
    final_survival_prob_excited_state = final_survival_prob[excited_state_index]
    print('final time:   ' + str(time_list[-1]) )

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    # ax.scatter( State_energy_for_ground_state , final_survival_prob_ground_state, label = 'electronic state |0>'  )
    # ax.scatter( State_energy_for_excited_state, final_survival_prob_excited_state, label = 'electronic state |1>' )

    ax.scatter(state_energy, final_survival_prob)

    ax.set_xlabel('energy')
    ax.set_ylabel('dilution factor')
    ax.set_title('dilution factor for all states')
    ax.set_yscale('log')
    # ax.legend(loc = 'best')

    # Plot how IPR change with time to see if IPR saturate or not
    fig1 = plt.figure(figsize = (10,10))
    spec1 = gridspec.GridSpec(nrows = 1, ncols = 1, figure = fig1)
    ax1 = fig1.add_subplot(spec1[0,0])

    state_energy_index = np.argsort(-state_energy )
    max_state_energy = state_energy[state_energy_index[0]]
    survival_prob_list_trans = np.transpose(survival_prob_list)

    for i in range(10):
        index = i
        if mode_number_list[state_energy_index[index]][0] == 0:
            label = '$|n_{e}> = $|0>'
        else:
            label = '$|n_{e}> = $|1>'
        label = label + " $|n_{v}>$ = " + str(mode_number_list[state_energy_index[index]][1:])
        label = label + " E = " + str(round( state_energy[state_energy_index[index]], 2))
        ax1.plot(time_list , survival_prob_list_trans[state_energy_index[index]] , label = label , linewidth = 2)
    ax1.legend(loc = 'best')

    ax1.set_xlabel('Time')
    ax1.set_ylabel('survival prob')
    ax1.set_yscale('log')
    ax1.set_title('survival prob for single state')


    if save_bool:
        fig_name = "energy vs dilution factor.png"
        fig_path = os.path.join(file_path , fig_name)
        fig.savefig(fig_path)

        fig1_name = "P(t) vs time.png"
        fig1_path = os.path.join(file_path, fig1_name)
        fig1.savefig(fig1_path)


def compare_survival_prob_subroutine(file_path_list, label_list, mode_index,  fig_file_path, save_bool):
    file_num = len(file_path_list)
    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    final_time = 1

    mode_number_list = []
    for i in range(file_num):
        file_path = file_path_list[i]
        label = label_list[i]
        state_num, nmode, state_energy, mode_number_list, survival_prob_list, time_list = Read_survival_prob_all_state(
            file_path)

        dt = time_list[1] - time_list[0]
        time_final_index = int(final_time / dt)

        ax.plot(time_list[:time_final_index], np.transpose(survival_prob_list)[mode_index][:time_final_index] , label = label , linewidth = 3)

    ax.legend(loc = 'best')
    ax.set_title('P(t) mode: ' + str(mode_number_list[mode_index][1:] ))
    ax.set_yscale('log')

    if save_bool:
        fig_name = "P(t) vs time diff file.png"
        fig_path = os.path.join(fig_file_path, fig_name)
        fig.savefig(fig_path)

def compute_dilution_factor(file_path):
    '''
    compute dilution factor from survival probability for all state.
    :return:
    '''
    state_num, nmode, state_energy, mode_number_list, survival_prob_list, time_list = Read_survival_prob_all_state(file_path)

    # average over second half time
    time_len_half = int(len(time_list)/2)

    survival_prob_list_trans = np.transpose(survival_prob_list)
    # row: state_num, col: time.
    dilution_factor = np.mean(survival_prob_list_trans[:,time_len_half:] , 1)

    return  state_energy, mode_number_list,dilution_factor
