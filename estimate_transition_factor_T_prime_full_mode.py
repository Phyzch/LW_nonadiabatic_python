'''
estimate nonadiabatic transition factor T' for 50 modes BChl. as function of Vt
'''
from util import *
from scipy.optimize import root
from Overlap_of_displaced_state import effective_num_coupling_submodule
from estimate_transition_factor_T_full_mode import  estimate_transition_energy,estimate_transition_factor_with_freq_and_energy_dimer, \
estimate_transition_factor_with_freq_and_temperature_dimer

def estimate_T_prime_prefactor_subroutine_dimer(EV_coupling_alpha_list, frequency_list, energy = 0, equal_energy_assumption = True):
    '''
     T_prime = Vt * prefactor.
     Here pre_factor is computed here.
    :param EV_coupling_alpha_list:
    :param frequency_list:
    :return:
    '''
    # Estimate local density of state Dq
    omega_rms = np.sqrt(np.mean(np.power(frequency_list, 2)))
    Dq = 1 / (np.pi * omega_rms)

    # compute Kt: connectivity
    dof = len(frequency_list)

    # at very low energy, we focus on Ki for quantum number n = 0.
    # energy here is energy for dimer.
    if equal_energy_assumption:
        # assume equal energy in each monomer
        qn = round( (energy/2) / (dof * omega_rms) )

        Ki_list = np.zeros([dof])
        for i in range(dof):
            alpha = EV_coupling_alpha_list[i]
            Ki = effective_num_coupling_submodule(alpha, qn)
            Ki_list[i] = Ki

        Kt = np.prod(Ki_list)

        Kt_dimer = np.power(Kt, 2)

    else:
        # assume energy concentrated in one monomer.
        qn1 = round(energy / (dof * omega_rms))

        Ki_list1 = np.zeros([dof])
        for i in range(dof):
            alpha = EV_coupling_alpha_list[i]
            Ki = effective_num_coupling_submodule(alpha, qn1)
            Ki_list1[i] = Ki
        Kt1 = np.prod(Ki_list1)

        Ki_list2 = np.zeros([dof])
        qn2 = 0
        for i in range(dof):
            alpha = EV_coupling_alpha_list[i]
            Ki = effective_num_coupling_submodule(alpha, qn2)
            Ki_list2[i] = Ki
        Kt2 = np.prod(Ki_list2)
        Kt_dimer = Kt1 * Kt2

    # estimate the nonadiabatic transition factor T_prime
    T_prime_prefactor = np.sqrt(2 * np.pi /3) * Dq * np.sqrt(Kt_dimer)

    return T_prime_prefactor

def estimate_T_prime_prefactor_temperature_subroutine_dimer(EV_coupling_alpha_list, frequency_list, temperature):
    '''

    :param EV_coupling_alpha_list:
    :param frequency_list:
    :param temperature:
    :return:
    '''
    # Estimate local density of state Dq
    omega_rms = np.sqrt(np.mean(np.power(frequency_list, 2)))
    Dq = 1 / (np.pi * omega_rms)

    # compute Kt: connectivity
    dof = len(frequency_list)

    # compute Kt
    quantum_number_list = np.array(range(10))
    qn_len = len(quantum_number_list)
    Ki_sqrt_list = []
    for mode_index in range(dof):
        frequency = frequency_list[mode_index]
        alpha = EV_coupling_alpha_list[mode_index]
        thermal_qn_prob = np.exp(-quantum_number_list * frequency/temperature) * (1 - np.exp(-frequency / temperature))
        Ki_n_list = []
        for j in range(qn_len):
            qn = quantum_number_list[j]
            Ki_n = effective_num_coupling_submodule(alpha,qn)
            Ki_n_list.append(Ki_n)

        Ki_n_sqrt_list = np.sqrt(Ki_n_list)
        Ki_sqrt = np.sum( thermal_qn_prob * Ki_n_sqrt_list ) # average over sqrt(Ki_n)
        Ki_sqrt_list.append(Ki_sqrt)

    Kn_sqrt = np.prod(Ki_sqrt_list)
    Kn_sqrt_dimer = Kn_sqrt * Kn_sqrt

    # estimate the nonadiabatic transition factor T_prime
    T_prime_prefactor = np.sqrt(2 * np.pi /3) * Dq * Kn_sqrt_dimer

    return T_prime_prefactor

def estimate_T_prime_prefactor_full_mode_dimer():
    '''

    :return:
    '''
    frequency_list = np.array([  84,  167,  183,  191,  214,  239,  256,  345,  368,  388,  407,
        423,  442,  473,  506,  565,  587,  623,  684,  696,  710,  727,
        776,  803,  845,  858,  890,  915,  967,  980, 1001, 1019, 1066,
       1089, 1105, 1117, 1137, 1158, 1180, 1190, 1211, 1229, 1252, 1289,
       1378, 1466, 1519, 1539, 1648, 1680])

    # data source : Table II of J. Chem. Phys. 134, 024506 (2011)
    Huang_Rhys_factor = 1/1000 * np.array([15.1, 8.1, 7.2, 19.6, 4.6, 7.8, 5.5, 16.1, 6.0, 4.1, 5.2, 2.9, 2.3, 1.7,
                                           1.6, 8.1, 3.9, 5.8, 2.3, 2.5, 1.8, 26.6, 9.5, 4.2, 2.5, 2.5,
                                           28.4, 4.8, 2.7, 3.1, 4.0, 9.7, 2.5, 2.1, 2.1, 10.3, 4.2,
                                           10.3, 2.5, 3.1, 2.5, 2.1, 2.5, 8.7, 8.6, 2.1, 1.7, 2.3, 2.5, 2.7])

    EV_coupling_alpha = np.sqrt(Huang_Rhys_factor)

    T_prime_prefactor_ground_state = estimate_T_prime_prefactor_subroutine_dimer(EV_coupling_alpha, frequency_list)

    return T_prime_prefactor_ground_state


def estimate_T_prime_prefactor_12_modes_largest_EV_dimer():
    '''

    :return:
    '''
    frequency_list = np.array([ 187, 335, 566, 730, 760, 895, 1020, 1115, 1131, 1162, 1281, 1376 ])

    # data source : Table II , III of J. Chem. Phys. 134, 024506 (2011)
    Huang_Rhys_factor = 1 / 1000 * np.array([19.6, 16.1, 8.1, 26.6, 9.5, 28.4, 9.7, 10.3, 4.2, 10.3, 8.7, 8.6])

    EV_coupling_alpha = np.sqrt(Huang_Rhys_factor)

    # assume energy E = 0.
    T_prime_prefactor_ground_state = estimate_T_prime_prefactor_subroutine_dimer(EV_coupling_alpha, frequency_list)

    return T_prime_prefactor_ground_state



def plot_E_scaling_factor_phase_diagram_for_BChl_dimer(frequency_list, Huang_Rhys_factor, scaling_factor_estimate, folder_path, save_bool):
    '''

    :return:
    '''
    matplotlib.rcParams.update({'font.size': 20})

    EV_coupling_alpha = np.sqrt(Huang_Rhys_factor)

    T_prime_prefactor = estimate_T_prime_prefactor_subroutine_dimer(EV_coupling_alpha, frequency_list)

    scaling_num = 20
    scaling_factor_list = np.linspace(0.08 , 0.3, scaling_num)

    # find transition energy when T_prime = 0
    T_prime = 0
    transition_energy_list_no_Vt = np.zeros([scaling_num])

    for i in range(scaling_num):
        scaling_factor = scaling_factor_list[i]
        result = root(estimate_transition_energy, x0 = np.array([1000]), args = (frequency_list, T_prime, scaling_factor) )
        transition_energy_no_Vt = result.x[0]
        transition_energy_list_no_Vt[i] = transition_energy_no_Vt

    # find transition energy when we turn Vt on
    Vt = 363
    T_prime = Vt * T_prime_prefactor
    print("T_prime  " + str(T_prime) + " for Vt: " + str(Vt))

    transition_energy_list_with_Vt = np.zeros([scaling_num])

    for i in range(scaling_num):
        scaling_factor = scaling_factor_list[i]
        result = root(estimate_transition_energy, x0 = np.array([1000]), args = (frequency_list, T_prime, scaling_factor))
        transition_energy_with_Vt = result.x[0]
        transition_energy_list_with_Vt[i] = transition_energy_with_Vt

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    ax.plot( 1/scaling_factor_list, transition_energy_list_no_Vt, marker = 'o' , linewidth = 2, label = 'Vt=0' )
    ax.plot(1/scaling_factor_list, transition_energy_list_with_Vt, marker = 'o' , linewidth = 2, label = 'Vt= 353 $cm^{-1}$')

    # plot the estimate for scaling factor using 1998 PNAS paper
    ax.axvline(1/scaling_factor_estimate, color = 'black', linewidth = 2)

    ax.legend(loc = 'best')
    ax.set_xlabel('1/$a$')
    ax.set_ylabel('Energy ($cm^{-1}$)')

    if save_bool:
        file_name = "transition_energy_E_vs_scaling_factor_dimer.svg"
        file_name = os.path.join(folder_path,file_name)
        fig.savefig(file_name)


    plt.show()


def plot_E_scaling_factor_phase_diagram_for_full_mode_BChl():
    '''

    :return:
    '''
    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/result 2022.10.06/paper fig/full_mode_T_T'_estimate/dimer_case/"
    frequency_list = np.array([  84,  167,  183,  191,  214,  239,  256,  345,  368,  388,  407,
        423,  442,  473,  506,  565,  587,  623,  684,  696,  710,  727,
        776,  803,  845,  858,  890,  915,  967,  980, 1001, 1019, 1066,
       1089, 1105, 1117, 1137, 1158, 1180, 1190, 1211, 1229, 1252, 1289,
       1378, 1466, 1519, 1539, 1648, 1680])

    # data source : Table II of J. Chem. Phys. 134, 024506 (2011)
    Huang_Rhys_factor = 1/1000 * np.array([15.1, 8.1, 7.2, 19.6, 4.6, 7.8, 5.5, 16.1, 6.0, 4.1, 5.2, 2.9, 2.3, 1.7,
                                           1.6, 8.1, 3.9, 5.8, 2.3, 2.5, 1.8, 26.6, 9.5, 4.2, 2.5, 2.5,
                                           28.4, 4.8, 2.7, 3.1, 4.0, 9.7, 2.5, 2.1, 2.1, 10.3, 4.2,
                                           10.3, 2.5, 3.1, 2.5, 2.1, 2.5, 8.7, 8.6, 2.1, 1.7, 2.3, 2.5, 2.7])
    dof = len(frequency_list)
    geometric_mean_freq = np.prod(np.power(frequency_list , 1/dof))
    scaling_factor_estimate = 1/270 * np.sqrt( geometric_mean_freq )

    plot_E_scaling_factor_phase_diagram_for_BChl_dimer(frequency_list, Huang_Rhys_factor, scaling_factor_estimate, folder_path, save_bool)

def plot_E_scaling_factor_phase_diagram_for_12_mode_BChl():
    '''

    :return:
    '''
    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/result 2022.10.06/paper fig/full_mode_T_T'_estimate/dimer_case/"
    frequency_list = np.array([ 187, 335, 566, 730, 760, 895, 1020, 1115, 1131, 1162, 1281, 1376 ])

    # data source : Table II , III of J. Chem. Phys. 134, 024506 (2011)
    Huang_Rhys_factor = 1/1000 *  np.array([19.6, 16.1, 8.1, 26.6, 9.5, 28.4, 9.7, 10.3, 4.2, 10.3, 8.7, 8.6])

    # manually input scaling factor list, estimate from mode type. for 12 mode model.
    scaling_factor_list = np.array([0.1, 0.14, 0.1, 0.1, 0.1, 0.13, 0.1, 0.15, 0.2, 0.17, 0.15, 0.15])
    dof = len(scaling_factor_list)
    scaling_factor_estimate = np.prod( np.power(scaling_factor_list, 1/dof) ) # geometric mean


    plot_E_scaling_factor_phase_diagram_for_BChl_dimer(frequency_list, Huang_Rhys_factor, scaling_factor_estimate, folder_path, save_bool)


def plot_E_transition_factor_subroutine(scaling_factor, Vt, Huang_Rhys_factor, frequency_list, fig_path, save_bool):
    '''

    :param V0: 3050 according to fitting.
    :param scaling_factor: scaling factor in Bigwood 1998 paper
    :param Vt:  nonadiabatic coupling strength
    :param Huang_Rhys_factor: EV coupling strength
    :return:
    '''
    matplotlib.rcParams.update({'font.size': 20})

    data_num = 20
    energy_list = np.linspace(0, 3000, data_num)
    Tq_list = np.zeros([data_num])
    T_prime_list =  np.zeros([data_num])

    EV_coupling_alpha = np.sqrt(Huang_Rhys_factor)


    for i in range(data_num):
        energy = energy_list[i]
        Tq = estimate_transition_factor_with_freq_and_energy_dimer(energy, frequency_list, scaling_factor)
        Tq_list[i] = Tq

        T_prime_prefactor = estimate_T_prime_prefactor_subroutine_dimer(EV_coupling_alpha, frequency_list, energy)
        T_prime = Vt * T_prime_prefactor
        T_prime_list[i] = T_prime

    Tq_T_prime_sum_list = Tq_list + T_prime_list

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])
    ax.plot(energy_list, Tq_list, linewidth=4)
    ax.plot(energy_list, Tq_T_prime_sum_list , linewidth = 4)

    ax.set_xlabel('E')
    ax.set_ylabel("T or (T + T')")
    # ax.axhline(y=1 , linewidth = 3, color = 'black')

    if save_bool:
        fig_name = "transition factor T, T' vs energy.svg"
        fig_name = os.path.join(fig_path, fig_name)
        fig.savefig(fig_name)

    plt.show()


def plot_temperature_transition_factor_subroutine(scaling_factor, Vt, Huang_Rhys_factor, frequency_list, fig_path, save_bool):
    '''

    :param V0: 3050 according to fitting.
    :param scaling_factor: scaling factor in Bigwood 1998 paper
    :param Vt:  nonadiabatic coupling strength
    :param Huang_Rhys_factor: EV coupling strength
    :return:
    '''
    matplotlib.rcParams.update({'font.size': 20})

    data_num = 20
    temperature_list = np.linspace(1, 3000, data_num)
    unit_transform = 208 / 300 # from K to cm^{-1}
    temperature_list_unit_wavenumber = unit_transform * temperature_list

    Tq_list = np.zeros([data_num])
    T_prime_list =  np.zeros([data_num])

    EV_coupling_alpha = np.sqrt(Huang_Rhys_factor)


    for i in range(data_num):
        temperature = temperature_list_unit_wavenumber[i]
        Tq = estimate_transition_factor_with_freq_and_temperature_dimer(temperature,frequency_list,scaling_factor)
        Tq_list[i] = Tq

        T_prime_prefactor = estimate_T_prime_prefactor_temperature_subroutine_dimer(EV_coupling_alpha, frequency_list, temperature)
        T_prime = Vt * T_prime_prefactor
        T_prime_list[i] = T_prime

    Tq_T_prime_sum_list = Tq_list + T_prime_list

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])
    ax.plot(temperature_list, Tq_list, linewidth=4)
    ax.plot(temperature_list, Tq_T_prime_sum_list , linewidth = 4)

    ax.set_xlabel('temperature(K)')
    ax.set_ylabel("T or (T + T')")

    # ax.axhline(y=1 , linewidth = 3, color = 'black')

    if save_bool:
        fig_name = "transition factor T, T' vs temperature.svg"
        fig_name = os.path.join(fig_path, fig_name)
        fig.savefig(fig_name)

    plt.show()

def plot_transition_factor_for_5_mode_BChl():
    '''

    :return:
    '''
    frequency_list = np.array([ 890, 727, 345, 1117, 1158 ])

    # data source : Table II , III of J. Chem. Phys. 134, 024506 (2011)
    Huang_Rhys_factor = 1/1000 *  np.array([28, 26.6, 16.1, 10.3, 10.3])

    # manually input scaling factor list, estimate from mode type. for 12 mode model.
    scaling_factor_list = np.array([0.110, 0.100, 0.069, 0.124, 0.126])
    dof = len(scaling_factor_list)
    scaling_factor_estimate = np.prod( np.power(scaling_factor_list, 1/dof) ) # geometric mean

    V0 = 3050
    Vt = 363

    fig_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/result 2022.10.06/paper fig/full_mode_T_T'_estimate/dimer_case_Bigwood formula/"
    save_bool = True
    # plot_E_transition_factor_subroutine(scaling_factor_estimate, Vt, Huang_Rhys_factor, frequency_list, fig_path, save_bool)

    plot_temperature_transition_factor_subroutine(scaling_factor_estimate, Vt, Huang_Rhys_factor, frequency_list, fig_path, save_bool)
def plot_transition_factor_for_full_mode_BChl():
    '''

    :return:
    '''
    save_bool = True
    fig_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/result 2022.10.06/paper fig/full_mode_T_T'_estimate/dimer_case_Bigwood formula/"
    frequency_list = np.array([  84,  167,  183,  191,  214,  239,  256,  345,  368,  388,  407,
        423,  442,  473,  506,  565,  587,  623,  684,  696,  710,  727,
        776,  803,  845,  858,  890,  915,  967,  980, 1001, 1019, 1066,
       1089, 1105, 1117, 1137, 1158, 1180, 1190, 1211, 1229, 1252, 1289,
       1378, 1466, 1519, 1539, 1648, 1680])

    # data source : Table II of J. Chem. Phys. 134, 024506 (2011)
    Huang_Rhys_factor = 1/1000 * np.array([15.1, 8.1, 7.2, 19.6, 4.6, 7.8, 5.5, 16.1, 6.0, 4.1, 5.2, 2.9, 2.3, 1.7,
                                           1.6, 8.1, 3.9, 5.8, 2.3, 2.5, 1.8, 26.6, 9.5, 4.2, 2.5, 2.5,
                                           28.4, 4.8, 2.7, 3.1, 4.0, 9.7, 2.5, 2.1, 2.1, 10.3, 4.2,
                                           10.3, 2.5, 3.1, 2.5, 2.1, 2.5, 8.7, 8.6, 2.1, 1.7, 2.3, 2.5, 2.7])
    dof = len(frequency_list)
    geometric_mean_freq = np.prod(np.power(frequency_list , 1/dof))
    scaling_factor_estimate = 1/270 * np.sqrt( geometric_mean_freq )

    Vt = 363

    # plot_E_transition_factor_subroutine(scaling_factor_estimate, Vt, Huang_Rhys_factor, frequency_list, fig_path, save_bool)

    plot_temperature_transition_factor_subroutine(scaling_factor_estimate, Vt, Huang_Rhys_factor, frequency_list, fig_path,
                                                  save_bool)

def compute_thermal_energy_BChl_full_mode():
    '''

    :return:
    '''
    frequency_list = np.array([  84,  167,  183,  191,  214,  239,  256,  345,  368,  388,  407,
        423,  442,  473,  506,  565,  587,  623,  684,  696,  710,  727,
        776,  803,  845,  858,  890,  915,  967,  980, 1001, 1019, 1066,
       1089, 1105, 1117, 1137, 1158, 1180, 1190, 1211, 1229, 1252, 1289,
       1378, 1466, 1519, 1539, 1648, 1680])

    T = 208.37 # 300K, in unit of cm^{-1}
    # T = 2500

    energy_in_each_mode = frequency_list / (np.exp(frequency_list / T) - 1)
    energy_in_monomer = np.sum(energy_in_each_mode)

    print("energy in monomer (full mode):" + str(energy_in_monomer))

    frequency_list = np.array([890, 727, 345, 1117, 1158])
    energy_in_each_mode = frequency_list / (np.exp(frequency_list / T) - 1)
    energy_in_monomer = np.sum(energy_in_each_mode)
    print("energy in monomer(5 mode): " + str(energy_in_monomer))



# plot_transition_factor_for_full_mode_BChl()
# plot_transition_factor_for_5_mode_BChl()
# compute_thermal_energy_BChl_full_mode()