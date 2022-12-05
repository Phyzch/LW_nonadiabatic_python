from util import *

def estimate_transition_factor_subroutine(Q, V0, scaling_factor,  frequency_list, energy ):
    '''
    V_{Q} = scaling_factor^{Q} * V0
    :param Q:
    :param V0:
    :param scaling_factor:
    :param energy:
    :return:
    '''
    dof = len(frequency_list)
    # connectivity K
    K = np.power(2 * dof, Q) /np.math.factorial(Q)

    # local density of state
    omega_rms = np.sqrt(np.mean(np.power(frequency_list, 2)))
    Dq = 1/ (omega_rms * np.power(Q,1/2) * np.pi)

    # coupling strength
    M = energy / (dof * omega_rms)
    Vq = V0 * np.power(scaling_factor , Q) * np.power(M, Q/2)

    Tq = np.sqrt(2 * np.pi / 3 ) * K * Dq * Vq

    return Tq

def estimate_transition_factor_with_freq_and_energy_dimer(energy, frequency_list, scaling_factor):
    '''
    choose formula for coupling Vq using Bigwood , Leitner, Gruebele, Wolynes 1998
        # here we assume equal energy in each monomer.

    :return:
    '''
    dof = len(frequency_list)

    # scaling factor: (use formula in 1998 PNAS Bigwood , Leitner, Gruebele, Wolynes)
    geometric_mean_freq = np.prod(np.power(frequency_list , 1/dof))
    # scaling_factor = 1/270 * np.sqrt( geometric_mean_freq )

    V0 = 3050
    # or specify V3
    # V3 = 2
    # V0 = V3 / np.power(scaling_factor ,3)

    #  here we assume equal energy in each monomer.
    monomer_energy = energy / 2

    # cubic order
    Q = 3
    Tq3 = estimate_transition_factor_subroutine(Q, V0, scaling_factor, frequency_list, monomer_energy )

    # quartic coupling
    Q = 4
    Tq4 = estimate_transition_factor_subroutine(Q, V0, scaling_factor, frequency_list, monomer_energy )

    Tq = Tq3 + Tq4

    dimer_Tq = Tq * 2

    return dimer_Tq

def estimate_transition_energy(energy, frequency_list , T_prime, scaling_factor):
    '''
    energy for Tq = 1
    :param energy:
    :param frequency_list:
    :return:
    '''
    Tq = estimate_transition_factor_with_freq_and_energy_dimer(energy, frequency_list, scaling_factor)

    return Tq + T_prime - 1

def estimate_transition_factor_full_mode_dimer():
    '''

    :return:
    '''
    matplotlib.rcParams.update({'font.size': 20})
    frequency_list = np.array([  84,  167,  183,  191,  214,  239,  256,  345,  368,  388,  407,
        423,  442,  473,  506,  565,  587,  623,  684,  696,  710,  727,
        776,  803,  845,  858,  890,  915,  967,  980, 1001, 1019, 1066,
       1089, 1105, 1117, 1137, 1158, 1180, 1190, 1211, 1229, 1252, 1289,
       1378, 1466, 1519, 1539, 1648, 1680])

    dof = len(frequency_list)
    geometric_mean_freq = np.prod(np.power(frequency_list , 1/dof))
    scaling_factor = 1/270 * np.sqrt( geometric_mean_freq )


    data_num = 20
    energy_list = np.linspace(0,2000, data_num)
    Tq_list = np.zeros([data_num])

    # here we assume equal energy in each monomer.

    for i in range(data_num):
        energy = energy_list[i]

        Tq = estimate_transition_factor_with_freq_and_energy_dimer(energy, frequency_list, scaling_factor)


        Tq_list[i] = Tq

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])
    ax.plot(energy_list , Tq_list, marker = 'o' , linewidth = 2)

    ax.set_xlabel('E')
    ax.set_ylabel('T')
    plt.show()

def estimate_transition_factor_12_most_strongly_coupled_mode():
    '''

    :return:
    '''
    matplotlib.rcParams.update({'font.size': 20})

    frequency_list = np.array([187, 335, 566, 730, 760, 895, 1020, 1115, 1131, 1162, 1281, 1376])

    scaling_factor_list = np.array([0.1, 0.14, 0.1, 0.1, 0.1, 0.13, 0.1, 0.15, 0.2, 0.17, 0.15, 0.15])
    dof = len(scaling_factor_list)

    scaling_factor = np.prod( np.power(scaling_factor_list, 1/dof) ) # geometric mean

    data_num = 20
    energy_list = np.linspace(0, 4000, data_num)
    Tq_list = np.zeros([data_num])

    for i in range(data_num):
        energy = energy_list[i]
        Tq = estimate_transition_factor_with_freq_and_energy_dimer(energy, frequency_list, scaling_factor)
        Tq_list[i] = Tq

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])
    ax.plot(energy_list, Tq_list, marker='o', linewidth=2)

    ax.set_xlabel('E')
    ax.set_ylabel('T')
    plt.show()
