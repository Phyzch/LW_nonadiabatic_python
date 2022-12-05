'''
estimate nonadiabatic transition factor T' for 50 modes BChl. as function of Vt
'''
from util import *
from scipy.optimize import root
from Overlap_of_displaced_state import effective_num_coupling_submodule
from estimate_transition_factor_T_full_mode import  estimate_transition_energy

def estimate_T_prime_prefactor_subroutine_dimer(EV_coupling_alpha_list, frequency_list):
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
    # at low energy, we focus on Ki for quantum number n = 0
    qn = 0
    Ki_list = np.zeros([dof])
    for i in range(dof):
        alpha = EV_coupling_alpha_list[i]
        Ki = effective_num_coupling_submodule(alpha, qn)
        Ki_list[i] = Ki

    Kt = np.prod(Ki_list)

    Kt_dimer = np.power(Kt, 2)

    # estimate the nonadiabatic transition factor T_prime
    T_prime_prefactor = np.sqrt(2 * np.pi /3) * Dq * np.sqrt(Kt_dimer)

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

    T_prime_prefactor = estimate_T_prime_prefactor_subroutine_dimer(EV_coupling_alpha, frequency_list)

    return T_prime_prefactor


def estimate_T_prime_prefactor_12_modes_largest_EV_dimer():
    '''

    :return:
    '''
    frequency_list = np.array([ 187, 335, 566, 730, 760, 895, 1020, 1115, 1131, 1162, 1281, 1376 ])

    # data source : Table II , III of J. Chem. Phys. 134, 024506 (2011)
    Huang_Rhys_factor = 1 / 1000 * np.array([19.6, 16.1, 8.1, 26.6, 9.5, 28.4, 9.7, 10.3, 4.2, 10.3, 8.7, 8.6])

    EV_coupling_alpha = np.sqrt(Huang_Rhys_factor)

    T_prime_prefactor = estimate_T_prime_prefactor_subroutine_dimer(EV_coupling_alpha, frequency_list)

    return T_prime_prefactor



def plot_E_scaling_factor_phase_diagram_for_BChl_dimer(frequency_list, Huang_Rhys_factor, scaling_factor_estimate):
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

    plt.show()


def plot_E_scaling_factor_phase_diagram_for_full_mode_BChl():
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
    dof = len(frequency_list)
    geometric_mean_freq = np.prod(np.power(frequency_list , 1/dof))
    scaling_factor_estimate = 1/270 * np.sqrt( geometric_mean_freq )

    plot_E_scaling_factor_phase_diagram_for_BChl_dimer(frequency_list, Huang_Rhys_factor, scaling_factor_estimate)

def plot_E_scaling_factor_phase_diagram_for_12_mode_BChl():
    '''

    :return:
    '''
    frequency_list = np.array([ 187, 335, 566, 730, 760, 895, 1020, 1115, 1131, 1162, 1281, 1376 ])

    # data source : Table II , III of J. Chem. Phys. 134, 024506 (2011)
    Huang_Rhys_factor = 1/1000 *  np.array([19.6, 16.1, 8.1, 26.6, 9.5, 28.4, 9.7, 10.3, 4.2, 10.3, 8.7, 8.6])

    # manually input scaling factor list, estimate from mode type. for 12 mode model.
    scaling_factor_list = np.array([0.1, 0.14, 0.1, 0.1, 0.1, 0.13, 0.1, 0.15, 0.2, 0.17, 0.15, 0.15])
    dof = len(scaling_factor_list)
    scaling_factor_estimate = np.prod( np.power(scaling_factor_list, 1/dof) ) # geometric mean


    plot_E_scaling_factor_phase_diagram_for_BChl_dimer(frequency_list, Huang_Rhys_factor, scaling_factor_estimate)

