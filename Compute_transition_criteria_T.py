import numpy as np

# This program is based on Leitner_Wolynes_1996: https://doi.org/10.1063/1.472920   eq. (4.6) to estimate T in molecules

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

frequency_list = [1149,508,291,474,843, 333]
quanta = [9, 2, 2, 1, 1, 1]

average_quanta = np.mean(quanta)

V0 = 300
scaling_factor = 0.15

T = T_in_same_electronic_state(frequency_list, V0, 1 / scaling_factor, average_quanta)

print(T)

