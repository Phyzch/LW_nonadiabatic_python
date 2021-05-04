import numpy as np
import matplotlib.pyplot as plt
# This program is based on Leitner_Wolynes_1996: https://doi.org/10.1063/1.472920   eq. (4.6) to estimate T in molecules
from Overlap_of_displaced_state import Overlap_n_m
def T_in_same_electronic_state(frequency, V0, inverse_scaling_factor, average_quanta):
    # N = dof
    N = len(frequency)

    frequency_rms =  np.sqrt(  np.sum ( np.power(frequency,2) ) / N )

    inverse_scaling_factor = inverse_scaling_factor / np.sqrt( max(average_quanta,1) )


    V3 = V0 * pow(1 / inverse_scaling_factor , 3)

    Phi_tilde = V3 * 1/ (np.pi * frequency_rms )

    kappa = 2 * N / inverse_scaling_factor

    F = 0
    for Q in range(3,5):
        F = F + pow (kappa / Q , Q) * np.exp(Q) /Q * np.power(2 * np.pi , -0.5)

    T = 2 * np.pi / 3 * pow(Phi_tilde,2) * pow(inverse_scaling_factor,6) * pow(F,2)

    return T

# frequency_list = [1149,508,291,474,843, 333]
# quanta = [9, 2, 2, 1, 1, 1]
# average_quanta = np.mean(quanta)
# V0 = 300
# scaling_factor = 0.15
# T = T_in_same_electronic_state(frequency_list, V0, 1 / scaling_factor, average_quanta)

def T_in_another_electronic_state(K, tunneling_strength, alpha, m, n, vibrational_frequency):
    # m is index for spin down. n is index for spin up
    D = 1 / vibrational_frequency
    Overlap = Overlap_n_m(alpha,m,n)
    coupling = tunneling_strength * Overlap

    T = 2 * np.pi / 3 * pow(K * coupling * D , 2)

    return T

def Analyze_T_scaling():
    frequency_list = [1149, 508, 291, 474, 843, 333]

    V0 = 300
    scaling_factor = [ 0.1 , 0.12, 0.15 , 0.18 , 0.2 ]

    scaling_num = len(scaling_factor)

    Crossing_point_state = [1, 0, 2, 2, 1, 1, 1 ]
    Crossing_point_complementary = [0, 9, 2, 2, 1, 1, 1 ]

    T_same_electronic_state = np.zeros([2, scaling_num])
    State_all = [Crossing_point_state, Crossing_point_complementary]
    for i in range(2):
        average_quanta = np.mean(State_all[i][1:])
        for j in range(scaling_num):
            scaling_factor_value = scaling_factor[j]
            T = T_in_same_electronic_state(frequency_list, V0, 1 / scaling_factor_value, average_quanta)
            T_same_electronic_state[i][j] = T

    color_list = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(scaling_factor, T_same_electronic_state[0], color=color_list[0], marker='o',
              label='spin up state:  ' + str(Crossing_point_state[1:]))
    ax.plot(scaling_factor, T_same_electronic_state[1], color=color_list[1], marker='o',
              label='spin down state:  ' + str(Crossing_point_complementary[1:]))
    ax.legend(loc='best')
    ax.set_xlabel(' scaling factor ')
    ax.set_ylabel('T in same spin state')
    ax.set_title('spin up state : ' + str(Crossing_point_state[1:]) + " \n spin down state:  " + str(
        Crossing_point_complementary[1:]))

    Criteria = [ (1 - T_same_electronic_state[0][i]) * (1 - T_same_electronic_state[1][i]) for i in range(scaling_num) ]
    fig1,ax1 = plt.subplots(nrows=1, ncols=1)
    ax1.plot(scaling_factor, Criteria , 'o')
    ax1.set_title('$(1-T_{1})$ * $(1 - T_{2})$')

    plt.show()

# Analyze_T_scaling()