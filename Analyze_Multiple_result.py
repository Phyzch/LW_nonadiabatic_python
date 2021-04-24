import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Generate_input_file import Generate_input_files
from Analyze_IPR import Analyze_IPR_all_state,Read_IPR_all_state, Read_local_density_of_state , Search_state


def Analyze_multiple_simulation_result():
    matplotlib.rcParams.update({'font.size': 15})
    parent_file_path = "/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/Batch_simulation_simple/"

    tunneling_strength = [0, 100 , 200 ]
    scaling_factor = [ 0.1 , 0.11, 0.12 , 0.13, 0.14]

    scaling_num = 5
    tunneling_coupling_num = 3

    # should in the form [ [same electronic_tunneling_strength] , [same electronic_tunneling_strength] , etc ]
    subfolder_path_list = []
    for i in range(tunneling_coupling_num):
        s1 = str(tunneling_strength[i])
        for j in range(scaling_num):
            s2 = str(scaling_factor[j])
            s3 = os.path.join(s1,s2)

            subfolder_path_list.append(s3)

    # subfolder_path_list = ['scaling_0.1_no_electronic_state_coupling', 'scaling_0.15_no_electronic_state_coupling', 'scaling_0.2_no_electronic_state_coupling',
    #                    'high_electronic_state_weak_coupling_0.1', 'high_electronic_state_weak_coupling_0.15', 'high_electronic_state_weak_coupling_0.2']
    #
    # scaling_num = 3
    # tunneling_coupling_num = 2

    file_num = len(subfolder_path_list)

    if(file_num != scaling_num * tunneling_coupling_num):
        raise NameError ("Error. file number not equal to scaling_num * tunneling_coupling_num ")

    path_list = []
    for i in range(file_num):
        path = os.path.join(parent_file_path, subfolder_path_list[i])
        path_list.append(path)

    shared_Mode_number_list = []
    shared_state_energy = []
    Final_IPR_list = []

    for i in range(file_num):
        state_num, nmode, state_energy, Mode_number_list, IPR_list, Time_list = Read_IPR_all_state(path_list[i])
        Final_IPR = IPR_list[-1]
        shared_Mode_number_list = Mode_number_list
        shared_state_energy = state_energy
        Final_IPR_list.append(Final_IPR)

    local_density_of_state_same_list = []
    local_density_of_state_another_list = []
    for i in range(file_num):
        state_num, local_density_of_state_same,\
        local_density_of_state_another, Mode_number_list = Read_local_density_of_state(path_list[i])

        local_density_of_state_same_list.append(local_density_of_state_same)
        local_density_of_state_another_list.append(local_density_of_state_another)


    Crossing_point_state = [1, 0, 2, 2, 1, 1, 1 ]
    Crossing_point_complementary = [0, 9, 2, 2, 1, 1, 1 ]

    position = Search_state(shared_Mode_number_list, Crossing_point_state)
    complementary_position = Search_state(shared_Mode_number_list, Crossing_point_complementary)
    if(position == -1 or complementary_position == -1):
        print("State or complementary state not in the list")
    else:
        # final IPR for state under different condition
        final_IPR_for_state = []
        complementary_final_IPR_for_state = []
        for i in range(file_num):
            final_IPR_for_state.append(Final_IPR_list[i][position])
            complementary_final_IPR_for_state.append(Final_IPR_list[i][complementary_position])

        # local density of state for specific state
        local_density_of_state_same_for_state = []
        local_density_of_state_another_for_state = []

        complementary_local_density_of_state_same_for_state =[]
        complementary_local_density_of_state_another_for_state = []

        for i in range(file_num):
            local_density_of_state_same_for_state.append(local_density_of_state_same_list[i][position])
            local_density_of_state_another_for_state.append(local_density_of_state_another_list[i][position])

            complementary_local_density_of_state_same_for_state.append(local_density_of_state_same_list[i][complementary_position])
            complementary_local_density_of_state_another_for_state.append(local_density_of_state_another_list[i][complementary_position])

        print(complementary_local_density_of_state_same_for_state)
        print(complementary_local_density_of_state_another_for_state)

        fig , ax = plt.subplots(nrows=1, ncols=1)
        color_list = ['blue' , 'red' ,'green' , 'purple']
        for i in range(tunneling_coupling_num):
            # with same tunneling, we should have same local_density_of_state_another in this case
            local_density_of_state_same_for_state_slice = local_density_of_state_same_for_state[i * scaling_num : (i+1) * scaling_num]
            common_local_density_of_state_another_for_state = round(
                local_density_of_state_another_for_state[i * scaling_num], 3)
            final_IPR_for_state_slice = final_IPR_for_state[i*scaling_num : (i+1) * scaling_num ]
            ax.plot( local_density_of_state_same_for_state_slice, final_IPR_for_state_slice, 'o', markersize = 10,
                     color=color_list[i],
                     label ='spin up  state: ' +str(Crossing_point_state[1:] ) + '  $N_{tunneling} = $' + str(common_local_density_of_state_another_for_state) )

            complementary_local_density_of_state_same_slice = \
                complementary_local_density_of_state_same_for_state[i * scaling_num : (i+1) * scaling_num]
            complementary_common_local_density_of_state_another_for_state = round(
                complementary_local_density_of_state_another_for_state[i * scaling_num], 3)

            complementary_final_IPR_for_state_slice = \
                complementary_final_IPR_for_state[i*scaling_num : (i+1) * scaling_num ]

            print(complementary_final_IPR_for_state_slice)

            # for clarity. we plot this against local density of state for same state.
            ax.plot(local_density_of_state_same_for_state_slice, complementary_final_IPR_for_state_slice, '*',  markersize = 10,
                    color = color_list[i],
                    label = 'spin down  state:' + str(Crossing_point_complementary[1:]) + '  $N_{tunneling} = $' + str(common_local_density_of_state_another_for_state))

        ax.legend(loc = 'best')
        ax.set_xlabel('$N_{anharmonicity}$   (spin up) ')
        ax.set_ylabel('IPR')
        ax.set_title('spin up state : ' +str(Crossing_point_state[1:]) + " \n spin down state:  "+ str(Crossing_point_complementary[1:]) )


        fig1, ax1 = plt.subplots(nrows=1,ncols=1)
        for i in range(tunneling_coupling_num):
            final_IPR_for_state_slice = final_IPR_for_state[i * scaling_num: (i + 1) * scaling_num]
            # T1
            T1 = np.array(local_density_of_state_same_for_state[
                                                          i * scaling_num: (i + 1) * scaling_num])
            # T1'
            T1_tilde =  np.array(local_density_of_state_another_for_state[i * scaling_num : (i + 1) * scaling_num ])

            # T2
            T2 = np.array (complementary_local_density_of_state_same_for_state[i * scaling_num: (i + 1) * scaling_num] )

            # T2'
            T2_tilde = np.array( complementary_local_density_of_state_another_for_state[i * scaling_num : (i + 1) * scaling_num ] )

            qualified_index = [ i for i in range(scaling_num) if T2[i] < 0.95 ]

            Criteria = T1_tilde[qualified_index] * T2_tilde[qualified_index] / ( (1-T1[qualified_index]) *  (1-T2[qualified_index]) )

            ax1.plot(Criteria, np.array(final_IPR_for_state_slice)[qualified_index] , 'o' , color = color_list[i] , label = 'spin up state  $ N_{tunneling} = $  ' + str(T1_tilde[0]))

        ax1.set_ylabel('IPR')
        ax1.set_xlabel("$ T'_{1} T'_{2} /[(1- T_{1}) * (1-T_{2})] $")
        ax1.legend(loc = 'best')

        plt.show()