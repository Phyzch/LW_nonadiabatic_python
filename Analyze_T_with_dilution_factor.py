from util import *
from Analyze_survival_prob import compute_dilution_factor
from anharmonic_transition_factor_T import estimate_anharmonic_transition_factor_T, estimate_anharmonic_transition_factor_T_method2
from nonadiabatic_transition_factor_T import analyze_nonadiabatic_transition_factor

def read_Vq(file_path):
    '''

    :param file_path:
    :return:
    '''
    file_name = os.path.join(file_path, "T_criterion.txt")
    with open(file_name) as f:
        data = f.read().splitlines()
        data = [i for i in data if i != '']

        line = re.split(' ', data[0])
        line = [int(i) for i in line]
        state_num, mode_num , energy_window = line

        mode_index_list = []
        Vq_list = []

        line_index = 0
        for i in range(int(state_num)):
            line_index = line_index + 1
            line = re.split(' ', data[line_index])
            mode_index = [int(i) for i in line if i != '']

            line_index = line_index + 1
            line = re.split(' ', data[line_index])
            line = [float(i) for i in line if i!='']
            _, Vq, _, _ = line

            mode_index_list.append(mode_index)
            Vq_list.append(Vq)

        mode_index_list = np.array(mode_index_list)
        Vq_list = np.array(Vq_list)

        return state_num, mode_num, mode_index_list, Vq_list

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

def estimate_anharmonic_transition_factor_T_state_lists_method2(Vq_list):
    T_list = []
    length = len(Vq_list)
    for i in range(length):
        Vq = Vq_list[i]
        T = estimate_anharmonic_transition_factor_T_method2(Vq)
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

def plot_dilution_factor_and_anharmonic_T():
    '''

    :return:
    '''
    # anharmonic coupling
    V0 = 300
    scaling_factor = 0.3

    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_phase_diagram/"
    path1 = "batch_simulation_dE=600/Vt=0,V0=300,a=0.3/"
    path2 = "batch_simulation_dE=600_high_energy/Vt=0/"
    path3 = "batch_simulation_dE=600_low_energy/Vt=0/"

    file_path_list = [path1, path2, path3]

    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    path_num = len(file_path_list)

    dilution_factor_list = []
    T_list = []
    mode_number_list_all_states = []
    for i in range(path_num):
        file_path = file_path_list[i]

        state_energy, mode_number_list, dilution_factor = compute_dilution_factor(file_path)

        dilution_factor = dilution_factor.tolist()

        mode_number_list = [ x[1:] for x in mode_number_list]

        # transition factor in same PES
        T = estimate_anharmonic_transition_factor_T_state_lists(V0, scaling_factor, mode_number_list)

        T_list = T_list + T
        dilution_factor_list = dilution_factor_list + dilution_factor
        mode_number_list_all_states = mode_number_list_all_states + mode_number_list

    T_list = np.array(T_list)
    dilution_factor_list = np.array(dilution_factor_list)
    mode_number_list_all_states = np.array(mode_number_list_all_states)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    ax.scatter(T_list, dilution_factor_list,  marker = 'o' , s=40 )

    ax.set_xlabel('T')
    ax.set_ylabel('$\sigma$')
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.show()

    if save_bool:
        fig_name = "anharmonic T vs dilution factor.svg"
        fig_name = os.path.join(folder_path, fig_name)
        fig.savefig(fig_name)

def analyze_dilution_factor_and_T_phsae_diagram():
    '''
    See section 2.3.2 in Leitner 2015 Quantum ergodicity and energy flow in molecules
    :return:
    '''
    # anharmonic coupling
    V0 = 300
    scaling_factor = 0.3

    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_phase_diagram/batch_simulation_dE=600/"

    # folder_path = '/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/study_ground_energy_shift_Vt=300/'
    file_path1 = 'Vt=0,V0=300,a=0.3'
    file_path2 = 'Vt=50,V0=300,a=0.3'
    file_path3 = 'Vt=100,V0=300,a=0.3'
    file_path4 = 'Vt=200,V0=300,a=0.3'
    file_path5 = 'Vt=300,V0=300,a=0.3'
    file_path6 = 'Vt=500,V0=300,a=0.3'

    file_path_list = [ file_path1, file_path2, file_path3, file_path4, file_path5, file_path6]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    # another folder contain states at higher energy
    folder_path1 = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_phase_diagram/batch_simulation_dE=600_high_energy/"
    file_path7 = "Vt=0"
    file_path8 = "Vt=100"
    file_path9 = "Vt=300"
    file_path10 = "Vt=500"
    file_path_list1 = [file_path7, file_path8, file_path9, file_path10]
    file_path_list1 = [os.path.join(folder_path1, path) for path in file_path_list1]

    file_path_list = file_path_list + file_path_list1

    folder_path2 = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_phase_diagram/batch_simulation_dE=600_low_energy/"
    file_path11 = "Vt=0"
    file_path12 = "Vt=100"
    file_path13 = "Vt=300"
    file_path14 = "Vt=500"
    file_path_list2 = [file_path11, file_path12, file_path13, file_path14]
    file_path_list2 = [os.path.join(folder_path2, path) for path in file_path_list2]
    file_path_list = file_path_list + file_path_list2

    Vt_list = [0, 50, 100, 200, 300, 500, 0 , 100, 300, 500,  0 ,100, 300, 500]

    path_num = len(file_path_list)

    T_list = []
    T_prime_list = []
    dilution_factor_list = []
    for i in range(path_num):
        file_path = file_path_list[i]
        Vt = Vt_list[i]

        state_energy, mode_number_list, dilution_factor = compute_dilution_factor(file_path)

        dilution_factor = dilution_factor.tolist()

        mode_number_list = [ x[1:] for x in mode_number_list]

        # transition factor in same PES
        T = estimate_anharmonic_transition_factor_T_state_lists(V0, scaling_factor, mode_number_list)

        # transition factor in different PES
        T_prime = analyze_nonadiabatic_transition_factor_lists(mode_number_list, Vt)

        # concatenate list
        T_list = T_list + T
        T_prime_list = T_prime_list + T_prime
        dilution_factor_list = dilution_factor_list + dilution_factor

    dilution_factor_list = np.array(dilution_factor_list)
    T_list = np.array(T_list)
    T_prime_list = np.array(T_prime_list)

    dilution_factor_criterion1 = 0.1
    dilution_factor_criterion2 = 0.01
    length = len(dilution_factor_list)

    localization_index_list = [index for index in range(length) if
                               dilution_factor_list[index] > dilution_factor_criterion1]
    mixing_index_list = [index for index in range(length) if
                         dilution_factor_list[index] <= dilution_factor_criterion1 and dilution_factor_list[index] > dilution_factor_criterion2 ]

    extended_state_index_list = [index for index in range(length) if
                                 dilution_factor_list[index] <= dilution_factor_criterion2]

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    ax.scatter(T_list[localization_index_list] , T_prime_list[localization_index_list] , marker = 'o' , color = 'red' , s=40)
    ax.scatter(T_list[extended_state_index_list] , T_prime_list[extended_state_index_list] , marker =  'v', color = 'blue', s=40)
    ax.scatter(T_list[mixing_index_list], T_prime_list[mixing_index_list], marker = 's' , color = 'purple', s = 40 )

    # sc = ax.scatter(T_list , T_prime_list , marker = 'o' ,  s=40 , c = np.log(dilution_factor_list), cmap = 'gray')
    # cbar = fig.colorbar(sc, ax = ax)
    # cbar.set_label("log( $\sigma$ )")

    T_boundary = np.linspace(0,1, 100)
    T_prime_boundary = 1 - T_boundary
    ax.plot(T_boundary, T_prime_boundary, color = 'black')

    ax.set_xlabel("$T$")
    ax.set_ylabel(" $T'$ ")
    ax.set_ylim([-0.05,0.5])

    if save_bool:
        fig_name = "dilution factor vs T.svg"
        fig_name = os.path.join(folder_path, fig_name)
        fig.savefig(fig_name)

    plt.show()