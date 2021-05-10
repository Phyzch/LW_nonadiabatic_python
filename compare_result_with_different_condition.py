import numpy as np
import os
import matplotlib.pyplot as plt
from Analyze_IPR import Analyze_IPR_all_state,Read_IPR_all_state, Read_local_density_of_state , Search_state, Read_average_coupling_strength_and_local_density_of_state
from Compute_transition_criteria_T import T_in_same_electronic_state, T_in_another_electronic_state, T_in_another_electronic_state_Duschrinsky_rotation
import matplotlib.gridspec as gridspec
import matplotlib

def Analyze_IPR_different_condition():
    matplotlib.rcParams.update({'font.size': 12})
    parent_file_path = "/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/"

    file_path_list = ["Duschrinsky_rotation/Batch_simulation_30/", "Duschrinsky_rotation/Batch_simulation_45/", "2_coordinate/Tunneling/Batch_simulation/" ,"1_coordinate/Tunneling/Batch_simulation_long_time/" ]

    Label_list = ["Duschrinsky rotation tunneling $\phi = 30$ " , "Duschrinsky rotation tunneling $\phi = 45$"  , "2D tunneling" , "1D tunneling"]

    frequency_list = [1149, 508, 291, 474, 843, 333]

    V0 = 300

    lambda_up = 2000
    lambda_down = 2000
    coupling_mode_freq = 1149
    alpha = (lambda_up + lambda_down) / coupling_mode_freq

    tunneling_strength = [ 0 , 10 , 20, 30, 50, 70 , 100, 200, 300, 500 ]
    scaling_factor = [  0.25 ]

    scaling_num = len(scaling_factor)
    tunneling_coupling_num = len(tunneling_strength)

    Crossing_point_state = [1, 3, 2, 2, 1, 1, 1]
    Crossing_point_complementary = [0, 3, 2, 2, 1, 1, 1]

    T_same_electronic_state = np.zeros([2, scaling_num])
    State_all = [Crossing_point_state, Crossing_point_complementary]
    for i in range(2):
        average_quanta = np.mean(State_all[i][1:])
        for j in range(scaling_num):
            scaling_factor_value = scaling_factor[j]
            T = T_in_same_electronic_state(frequency_list, V0, 1 / scaling_factor_value, average_quanta)
            T_same_electronic_state[i][j] = T

    color_list = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    marker_list = ['*', 'o', 's','p']

    fig = plt.figure(figsize=(10, 20))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    spec.update(hspace=0.5, wspace=0.3)
    ax_1 = fig.add_subplot(spec[0, 0])

    Len = len(file_path_list)
    # should in the form [ [same electronic_tunneling_strength] , [same electronic_tunneling_strength] , etc ]
    for condition_index in range(Len):
        Folder_path = os.path.join(parent_file_path, file_path_list[condition_index])
        subfolder_path_list = []
        for i in range(tunneling_coupling_num):
            s1 = str(tunneling_strength[i])
            for j in range(scaling_num):
                s2 = str(scaling_factor[j])
                s3 = os.path.join(s1,s2)

                subfolder_path_list.append(s3)

        file_num = len(subfolder_path_list)

        if(file_num != scaling_num * tunneling_coupling_num):
            raise NameError ("Error. file number not equal to scaling_num * tunneling_coupling_num ")

        path_list = []
        for i in range(file_num):
            path = os.path.join(Folder_path, subfolder_path_list[i])
            path_list.append(path)

        shared_Mode_number_list = []
        shared_state_energy = []
        Final_IPR_list = []

        for i in range(file_num):
            state_num, nmode, state_energy, Mode_number_list, IPR_list, Time_list = Read_IPR_all_state(path_list[i])
            Len = len(IPR_list)
            # dilution_factor = 1/np.array(IPR_list)
            # average_dilution_factor = np.mean(dilution_factor[int(Len/2) :] , 0)
            # Final_IPR = 1/ average_dilution_factor
            Final_IPR = np.mean(IPR_list[int(Len  / 5) : ] , 0)
            shared_Mode_number_list = Mode_number_list
            shared_state_energy = state_energy
            Final_IPR_list.append(Final_IPR)


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

            ax_1.plot(np.array(tunneling_strength) + 0.1, final_IPR_for_state, marker = marker_list[condition_index], markersize = 10, color = color_list[condition_index] , label = Label_list[condition_index])



    ax_1.set_xlabel('Tunneling strength t')

    ax_1.legend(loc = 'best')

    ax_1.set_ylabel('IPR')

    ax_1.set_xscale('log')

    ax_1.set_title('state:  ' + str(Crossing_point_state[1:])  + "  T =  " + str( round( T_same_electronic_state[0][0] , 3 ) ))

    plt.show()