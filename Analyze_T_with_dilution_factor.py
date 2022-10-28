from util import *
from Analyze_survival_prob import compute_dilution_factor
from anharmonic_transition_factor_T import estimate_anharmonic_transition_factor_T
from nonadiabatic_transition_factor_T import analyze_nonadiabatic_transition_factor

def estimate_anharmonic_transition_factor_T_state_lists(V0, scaling_factor, mode_number_list):
    '''

    :param V0:
    :param scaling_factor:
    :param mode_number_list:
    :return:
    '''
    T_list = []
    length = len(mode_number_list)
    for i in range(length):
        mode_number = mode_number_list[i]
        T = estimate_anharmonic_transition_factor_T(V0, scaling_factor, mode_number)
        T_list.append(T)

    return T_list

def analyze_nonadiabatic_transition_factor_lists(mode_number_list, Vt):
    '''

    :param mode_number_list:
    :param Vt:
    :return:
    '''
    T_prime_list = []
    length = len(mode_number_list)
    for i in range(length):
        mode_number = mode_number_list[i]
        T_prime = analyze_nonadiabatic_transition_factor(mode_number, Vt)
        T_prime_list.append(T_prime)

    return T_prime_list

def analyze_dilution_factor_and_T():
    '''

    :return:
    '''
    # anharmonic coupling
    V0 = 300
    scaling_factor = 0.3

    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_dE=600/"

    # folder_path = '/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/study_ground_energy_shift_Vt=300/'
    file_path1 = 'Vt=0,V0=300,a=0.3'
    file_path2 = 'Vt=50,V0=300,a=0.3'
    file_path3 = 'Vt=100,V0=300,a=0.3'
    file_path4 = 'Vt=200,V0=300,a=0.3'
    file_path5 = 'Vt=300,V0=300,a=0.3'
    file_path6 = 'Vt=500,V0=300,a=0.3'

    Vt_list = [0, 50, 100, 200, 300, 500]

    file_path_list = [ file_path1, file_path2, file_path3, file_path4, file_path5, file_path6]

    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    path_num = len(file_path_list)

    T_list = []
    T_prime_list = []
    dilution_factor_list = []
    for i in range(path_num):
        file_path = file_path_list[i]
        Vt = Vt_list[i]

        state_energy, mode_number_list, dilution_factor = compute_dilution_factor(file_path)

        mode_number_list = [ x[1:] for x in mode_number_list]

        # transition factor in same PES
        T = estimate_anharmonic_transition_factor_T_state_lists(V0, scaling_factor, mode_number_list)
        # transition factor in different PES
        T_prime = analyze_nonadiabatic_transition_factor_lists(mode_number_list, Vt)

        T_list.append(T)
        T_prime_list.append(T_prime)
        dilution_factor_list.append(dilution_factor)

    T_list_flat = np.array(T_list).flatten()
    T_prime_list_flat = np.array(T_prime_list).flatten()
    dilution_factor_list_flat = np.array(dilution_factor_list).flatten()

    dilution_factor_criterion = 0.05
    length = len(dilution_factor_list_flat)

    localization_index_list = [index  for index in range(length)  if dilution_factor_list_flat[index] > dilution_factor_criterion]
    extended_state_index_list = [index for index in range(length) if dilution_factor_list_flat[index] <= dilution_factor_criterion ]

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    ax.scatter(T_list_flat[localization_index_list] , T_prime_list_flat[localization_index_list] , marker = 'o' , color = 'red' , s=40)
    ax.scatter(T_list_flat[extended_state_index_list] , T_prime_list_flat[extended_state_index_list] , marker =  'v', color = 'blue', s=40)

    T_boundary = np.linspace(0,1, 100)
    T_prime_boundary = 1 - T_boundary
    ax.plot(T_boundary, T_prime_boundary, color = 'black')

    ax.set_xlabel("$T_{anharmonic}$")
    ax.set_ylabel(" $T_{nonadiabatic}$ ")
    plt.show()