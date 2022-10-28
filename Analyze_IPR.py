from util import *


def Read_IPR_all_state(file_path):
    file_name = os.path.join(file_path,"IPR_all_state.txt")

    state_num = 0
    state_energy = []
    Mode_number_list = []
    Time_list = []
    IPR_list = []
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

            Mode_number_list.append(mode_number)
            line_index = line_index + 1

        step_num = int ((datalen - line_index)/2 )

        for i in range(step_num):
            line = data[line_index]
            line = re.split(' ', line)

            t = float(line[0])
            Time_list.append(t)

            line_index = line_index  + 1

            line = data[line_index]
            line = re.split(' ', line)

            IPR = [float(i) for i in line if i!= '']

            IPR_list.append(IPR)

            line_index = line_index + 1

        state_energy = np.array(state_energy)
        Mode_number_list = np.array(Mode_number_list)
        IPR_list = np.array(IPR_list )
        Time_list = np.array(Time_list)

        nmode = len(Mode_number_list[0])
        return state_num, nmode, state_energy, Mode_number_list, IPR_list, Time_list

def Read_local_density_of_state(file_path):
    file_name = os.path.join(file_path, "local_density_of_state.txt" )

    with open(file_name) as f:
        data = f.read()
        data = re.split('\n', data)
        data = [i for i in data if i != '']
        datalen = len(data)

        line = data[0]
        line = re.split(' ', line)
        state_num = int(line[0])

        line = data[1]
        line = re.split(' ',line)
        local_density_of_state_same_electronic_state = [float(i) for i in line if i!='']

        line = data[2]
        line =re.split(' ',line)
        local_density_of_state_in_another_electronic_state = [float(i) for i in line if i!='']

        line_index = 3
        Mode_number_list = []
        for i in range(state_num):
            line = data[line_index]
            line = re.split(' ',line)
            mode_num = [int(i) for i in line if i!='']
            Mode_number_list.append(mode_num)
            line_index = line_index + 1

        return state_num, local_density_of_state_same_electronic_state, local_density_of_state_in_another_electronic_state, Mode_number_list

def Read_average_coupling_strength_and_local_density_of_state(folder_path):
    file_name = os.path.join(folder_path, "transition_criteria_info.txt")

    with open(file_name) as f:
        data = f.readlines()
        datalen = len(data)

        line = data[0].strip('\n')
        line = re.split(' ', line)
        state_num = int(line[0])

        line = data[1].strip('\n')
        line = re.split(' ',line)
        rho_local_same = [float(i) for i in line if i!='']
        rho_local_same = np.array(rho_local_same)

        line = data[2].strip('\n')
        line = re.split(' ',line)
        rho_local_another = [float(i) for i in line if i!='']
        rho_local_another = np.array(rho_local_another)

        line = data[3].strip('\n')
        line = re.split(' ',line)
        V_average_same = [float(i) for i in line if i!='']
        V_average_same = np.array(V_average_same)

        line = data[4].strip('\n')
        line = re.split(' ',line)
        V_average_another = [float(i) for i in line if i!='']
        V_average_another = np.array(V_average_another)

        line_index = 5
        Mode_number_list = []
        for i in range(state_num):
            line = data[line_index].strip('\n')
            line = re.split(' ',line)
            mode_num = [int(i) for i in line if i!='']
            Mode_number_list.append(mode_num)
            line_index = line_index + 1

        Criteria_T_same = 2*np.pi / 3 * pow(V_average_same * rho_local_same,2)
        Criteria_T_another = 2 * np.pi /3 * pow(V_average_another * rho_local_another ,2 )

        return Mode_number_list , rho_local_same, rho_local_another, V_average_same, V_average_another , Criteria_T_same, Criteria_T_another



def Search_state(mode_state_list, mode_state):
    Len = len(mode_state_list)
    mode_state_list_1 = mode_state_list.tolist()
    for i in range(Len):
        if(mode_state_list_1[i] == mode_state):
            return i

    return -1

def Analyze_IPR_all_state(file_path , save_bool):
    matplotlib.rcParams.update({'font.size': 15})
    state_num, nmode , state_energy, mode_number_list, IPR_list , time_list = Read_IPR_all_state(file_path)

    Ground_state_index = [ i for i in range(state_num) if mode_number_list[i][0] == 0 ]
    Excited_state_index = [i for i in range(state_num) if mode_number_list[i][0] == 1 ]

    State_energy_for_ground_state = state_energy[Ground_state_index]
    State_energy_for_excited_state = state_energy[Excited_state_index]

    # take IPR at final time as final IPR
    time_len = len(time_list)
    start_time_to_plot_index = int(time_len * 2 / 3)

    Final_IPR = np.mean(IPR_list[ start_time_to_plot_index: ] , 0)
    Final_IPR_for_ground_state = Final_IPR[Ground_state_index]
    Final_IPR_for_excited_state = Final_IPR[Excited_state_index]
    print('final time:   ' + str(time_list[-1]) )

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    # ax.scatter( State_energy_for_ground_state , Final_IPR_for_ground_state, label = 'electronic state |0>'  )
    # ax.scatter( State_energy_for_excited_state, Final_IPR_for_excited_state, label = 'electronic state |1>' )

    ax.scatter(state_energy, Final_IPR)

    ax.set_xlabel('energy')
    ax.set_ylabel('IPR')
    ax.set_title('IPR for all states')
    # ax.set_ylim([0,25])
    # ax.legend(loc = 'best')

    # Plot how IPR change with time to see if IPR saturate or not
    fig1 = plt.figure(figsize = (10,10))
    spec1 = gridspec.GridSpec(nrows = 1, ncols = 1, figure = fig1)
    ax1 = fig1.add_subplot(spec1[0,0])

    state_energy_index = np.argsort(-state_energy )
    max_state_energy = state_energy[state_energy_index[0]]
    IPR_list_trans = np.transpose(IPR_list)

    for i in range(1):
        index = i
        if mode_number_list[state_energy_index[index]][0] == 0:
            label = '$|n_{e}> = $|0>'
        else:
            label = '$|n_{e}> = $|1>'
        label = label + " $|n_{v}>$ = " + str(mode_number_list[state_energy_index[index]][1:])
        label = label + " E = " + str(round( state_energy[state_energy_index[index]], 2))
        ax1.plot(time_list , IPR_list_trans[state_energy_index[index]] , label = label , linewidth = 2)
    ax1.legend(loc = 'best')

    ax1.set_xlabel('Time')
    ax1.set_ylabel('IPR')
    ax1.set_title('IPR for single state')
    # ax1.set_ylim([0,60])


    if save_bool:
        fig_name = "energy vs IPR.png"
        fig_path = os.path.join(file_path , fig_name)
        fig.savefig(fig_path)

        fig1_name = "IPR vs time.png"
        fig1_path = os.path.join(file_path, fig1_name)
        fig1.savefig(fig1_path)

def compare_IPR_subroutine(file_path_list, label_list, mode_index,  fig_file_path, save_bool):
    file_num = len(file_path_list)
    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    mode_number_list = []

    for i in range(file_num):
        file_path = file_path_list[i]
        label = label_list[i]
        state_num, nmode, state_energy, mode_number_list, IPR_list, time_list = Read_IPR_all_state(
            file_path)

        ax.plot(time_list, np.transpose(IPR_list)[mode_index] , label = label , linewidth = 3)

    ax.legend(loc = 'best')
    ax.set_title('IPR mode: ' + str(mode_number_list[mode_index]))
    ax.set_yscale('log')

    if save_bool:
        fig_name = "IPR diff file.png"
        fig_path = os.path.join(fig_file_path, fig_name)
        fig.savefig(fig_path)
