import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Generate_input_file import Generate_input_files
from Analyze_IPR import Analyze_IPR_all_state,Read_IPR_all_state, Read_local_density_of_state , Search_state, Read_average_coupling_strength_and_local_density_of_state
from Compute_transition_criteria_T import T_in_same_electronic_state, T_in_another_electronic_state, T_in_another_electronic_state_Duschrinsky_rotation
import matplotlib.gridspec as gridspec


def Analyze_multiple_simulation_result():
    matplotlib.rcParams.update({'font.size': 15})
    parent_file_path = "/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/2_coordinate/Tunneling/Batch_simulation/"

    frequency_list = [1149, 508, 291, 474, 843, 333]

    V0 = 300

    lambda_up = 2000
    lambda_down = 2000
    coupling_mode_freq = 1149
    alpha = (lambda_up + lambda_down) / coupling_mode_freq

    tunneling_strength = [ 0 , 10 , 20, 30, 50, 70 , 100, 200, 300 , 500 ]
    scaling_factor = [ 0.1,  0.15, 0.16, 0.17, 0.18, 0.19, 0.2 ]

    scaling_num = len(scaling_factor)
    tunneling_coupling_num = len(tunneling_strength)

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
        Len = len(IPR_list)
        Final_IPR = np.mean(IPR_list[int(Len  / 2) : ] , 0)
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

    # for i in range(file_num):
    #     Mode_number_list, rho_local_same, rho_local_another, V_average_same, V_average_another, \
    #     Criteria_T_same, Criteria_T_another \
    #     = Read_average_coupling_strength_and_local_density_of_state(path_list[i])


    Crossing_point_state = [1, 3, 2, 2, 1, 1, 1 ]
    Crossing_point_complementary = [0, 3, 2, 2, 1, 1, 1 ]

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

        # compute T:
        T_same_electronic_state = np.zeros([2,scaling_num])
        State_all = [Crossing_point_state, Crossing_point_complementary]
        for i in range(2):
            average_quanta = np.mean(State_all[i][1:])
            for j in range (scaling_num):
                scaling_factor_value = scaling_factor[j]
                T = T_in_same_electronic_state(frequency_list, V0, 1/scaling_factor_value, average_quanta)
                T_same_electronic_state[i][j] = T

        # compute T':
        T_tunneling = np.zeros([2, tunneling_coupling_num])
        # connectivity for coupling K
        K  = 100
        state_m = Crossing_point_complementary[1]
        state_n = Crossing_point_state[1]
        for i in range(2):
            for j in range(tunneling_coupling_num):
                # T_prime = T_in_another_electronic_state(K,tunneling_strength[j],alpha, state_m, state_n, coupling_mode_freq )

                # For Duschrinsky rotation .
                typical_state_overlap = 0.08
                T_prime = T_in_another_electronic_state_Duschrinsky_rotation(K,tunneling_strength[j],typical_state_overlap, frequency_list)

                T_tunneling[i][j] = T_prime

        fig = plt.figure(figsize=(10, 20))
        spec = gridspec.GridSpec(nrows=1, ncols=2, figure=fig)
        spec.update(hspace=0.5, wspace=0.3)
        ax_1 = fig.add_subplot(spec[0,0])
        ax_2 = fig.add_subplot(spec[0,1])

        color_list = ['blue' , 'orange' ,'green' , 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
        for i in range(tunneling_coupling_num):
            # with same tunneling, we should have same local_density_of_state_another in this case

            final_IPR_for_state_slice = final_IPR_for_state[i*scaling_num : (i+1) * scaling_num ]
            ax_1.plot( scaling_factor, final_IPR_for_state_slice, 'o', markersize = 10,
                     color=color_list[i],
                     label ='spin up  state: ' +str(Crossing_point_state[1:] ) + '  $T_{tunneling} = $' + str(round(T_tunneling[0][i] , 4 )) )

            complementary_final_IPR_for_state_slice = complementary_final_IPR_for_state[ i * scaling_num : (i+1) * scaling_num ]

            # for clarity. we plot this against local density of state for same state.
            ax_1.plot( scaling_factor, complementary_final_IPR_for_state_slice, '*',  markersize = 10,
                    color = color_list[i],
                    label = 'spin down  state:' + str(Crossing_point_complementary[1:]) + '  $T_{tunneling} = $' + str( round(T_tunneling[0][i] , 4 )) )


        ax_1.legend(loc = 'best', prop={'size': 8})
        ax_1.set_xlabel(' scaling factor ')
        ax_1.set_ylabel('IPR')
        ax_1.set_title('spin up state : ' +str(Crossing_point_state[1:]) + " \n spin down state:  "+ str(Crossing_point_complementary[1:]) )

        ax_2.plot(scaling_factor, T_same_electronic_state[0], color = color_list[0], marker = 'o', label = 'spin up state:  '+str(Crossing_point_state[1:] ))
        ax_2.plot(scaling_factor, T_same_electronic_state[1], color = color_list[1], marker = 'o', label = 'spin down state:  ' + str(Crossing_point_complementary[1:]) )
        ax_2.legend(loc = 'best')
        ax_2.set_xlabel(' scaling factor ')
        ax_2.set_ylabel('T in same spin state')
        ax_2.set_title('spin up state : ' +str(Crossing_point_state[1:]) + " \n spin down state:  "+ str(Crossing_point_complementary[1:]))


        fig1, ax1 = plt.subplots(nrows=1,ncols=1)
        fig2, ax2 = plt.subplots(nrows=1, ncols=1)
        C_fitting = 1.1
        for i in range(scaling_num):
            final_IPR_for_state_slice = [final_IPR_for_state[i + j * scaling_num] for j in range(tunneling_coupling_num)]
            final_IPR_for_complementary_state_slice = complementary_final_IPR_for_state[i * scaling_num: (i + 1) * scaling_num]
            # T1
            T1 = T_same_electronic_state[0][i]
            # T1'
            T1_tilde =  T_tunneling[0]

            # T2
            T2 = T_same_electronic_state[1][i]

            # T2'
            T2_tilde = T_tunneling[1]


            if(T1 < C_fitting and T2 < C_fitting):
                qualified = True
            else:
                qualified = False

            if qualified:
                for j in range(tunneling_coupling_num):
                    if(T1_tilde[j] == 0 or T2_tilde[j] == 0):
                        T1_tilde[j] = pow(10,-4)
                        T2_tilde[j] = pow(10,-4)

                Criteria_1 = T1_tilde  * T2_tilde / ( (C_fitting-T1) *  (C_fitting-T2) )
                Criteria = T1_tilde * T2_tilde
                # plot different T' with different color
                # ax1.plot(Criteria, np.array(final_IPR_for_state_slice)[qualified_index] , 'o' , color = color_list[i] , label = '$T_{t,up}$ = ' + str( round(T1_tilde, 3))
                #          + "  $T_{t,down}$ =  " + str( round(T2_tilde,3)  ) )
                #
                # ax1.plot(Criteria, np.array(final_IPR_for_complementary_state_slice)[qualified_index], '*', markersize = 10, color=color_list[i])


                # plot different scaling with different color
                Criteria_len = len(Criteria)
                final_IPR_for_state_slice_qualified = np.array(final_IPR_for_state_slice)
                final_IPR_for_complementary_state_qualified = np.array(final_IPR_for_complementary_state_slice)

                ax1.plot(Criteria, final_IPR_for_state_slice_qualified, marker= 'o', color=color_list[i] , label = '$T_{up}$ = ' + str( round (T_same_electronic_state[0][i] , 3 ) )
                          + ' $T_{down}$ = ' + str(  round (T_same_electronic_state[1][i] , 3  )  ) )


                ax2.plot(Criteria_1 , final_IPR_for_state_slice_qualified, marker = 'o', color=color_list[i] , label = '$T_{up}$ = ' + str( round (T_same_electronic_state[0][i] , 3 ) )
                          + ' $T_{down}$ = ' + str(  round (T_same_electronic_state[1][i] , 3  )  ))



        ax1.set_ylabel('IPR')
        # ax1.set_xlabel("$ T'_{1} T'_{2} /[(C- T_{1}) * (C-T_{2})] $   C = " + str(C_fitting))
        ax1.set_xlabel("$ T'_{1} T'_{2}  $")
        ax1.set_xscale('log')
        ax1.set_xticks([pow(10,-7) , pow(10,-5), pow(10,-3) , pow(10,-1) , 1, pow(10,1) ])
        # ax1.set_yscale('log')
        ax1.legend(loc = 'best')


        ax2.set_ylabel('IPR')
        ax2.set_xlabel("$T'_{1}  T'_{2} / (1-T_{1}) (1- T_{2})$")
        ax2.legend(loc = 'best')
        ax2.set_xscale('log')
        ax2.set_xticks([pow(10,-7) , pow(10,-5), pow(10,-3) , pow(10,-1) , 1, pow(10,1) , pow(10,3) ])

        plt.show()