import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
import os

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

        Step_num = int ((datalen - line_index)/2 )

        for i in range(Step_num):
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

def Analyze_IPR_all_state(file_path):
    matplotlib.rcParams.update({'font.size': 15})
    state_num, nmode , state_energy, Mode_number_list, IPR_list , Time_list = Read_IPR_all_state(file_path)

    Ground_state_index = [ i for i in range(state_num) if Mode_number_list[i][0] == 0 ]
    Excited_state_index = [i for i in range(state_num) if Mode_number_list[i][0] == 1 ]

    State_energy_for_ground_state = state_energy[Ground_state_index]
    State_energy_for_excited_state = state_energy[Excited_state_index]

    Final_IPR = IPR_list[-1]
    Final_IPR_for_ground_state = Final_IPR[Ground_state_index]
    Final_IPR_for_excited_state = Final_IPR[Excited_state_index]
    print('final time:   ' + str(Time_list[-1]) )

    fig,ax = plt.subplots(nrows=1, ncols=1)
    ax.scatter( State_energy_for_ground_state , Final_IPR_for_ground_state, label = 'ground state' )
    ax.scatter( State_energy_for_excited_state, Final_IPR_for_excited_state, label = 'excited state')

    ax.set_xlabel('energy')
    ax.set_ylabel('IPR')
    ax.set_title('IPR for all states')
    ax.legend(loc = 'best')

    # Plot how IPR change with time to see if IPR saturate or not
    fig1, ax1 = plt.subplots(nrows=1, ncols=1)

    state_energy_index = np.argsort(-state_energy )
    max_state_energy = state_energy[state_energy_index[0]]
    IPR_list_trans = np.transpose(IPR_list)

    ax1.plot(Time_list , IPR_list_trans[state_energy_index[0]] )
    ax1.set_xlabel('Time')
    ax1.set_ylabel('IPR')
    ax1.set_title('IPR for single state')

    plt.show()