import numpy

from Overlap_of_displaced_state import effective_num_coupling_submodule
from util import *
from Analyze_Nloc_with_dilution_factor import read_Nloc_anharmonic_nonadiabatic_coupling

def analyze_nonadiabatic_transition_factor(qn_list, Vt):
    '''
    Vt: nonadiabatic coupling strength
    qn_list : quantum number
    :return:
    '''
    alpha_list = np.array([ 0.169, 0.163, 0.127, 0.101, 0.101])

    mode_num = len(qn_list)
    effective_num_list = np.zeros([mode_num])
    for i in range(mode_num):
        alpha = alpha_list[i]
        qn = qn_list[i]
        effective_num = effective_num_coupling_submodule(alpha, qn)
        effective_num_list[i] = effective_num

    # print("effective_num_each_mode : " + str(effective_num_list))
    effective_num_prod = np.prod(effective_num_list)
    # print("effective coupling num : " + str(effective_num_prod))

    # compute nonadiabatic transition factor T
    frequency_list = np.array([890, 727, 345, 1117, 1158])
    D = 1 / (np.pi * np.mean(frequency_list))
    K = effective_num_prod

    # Franck condon factor is approximated as 1/sqrt(K)
    T = np.sqrt(2 * np.pi / 3) * D * Vt * np.sqrt(K)

    # print("state : " + str(qn_list))
    # print("Vt = " + str(Vt) +  " transition criterion $T'$ : " + str(T))

    return T

def compare_nonadiabatic_T_and_Nloc():
    '''

    :return:
    '''
    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_dE=600/"

    file_path2 = 'Vt=50,V0=300,a=0.3'
    file_path3 = 'Vt=100,V0=300,a=0.3'
    file_path4 = 'Vt=200,V0=300,a=0.3'
    file_path5 = 'Vt=300,V0=300,a=0.3'
    file_path6 = 'Vt=500,V0=300,a=0.3'

    file_path_list = [ file_path2, file_path3, file_path4, file_path5, file_path6]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    Vt_list = [ 50, 100, 200, 300 , 500]

    path_num = len(file_path_list)

    diff_PES_Nloc_list = []
    nonadiabatic_T_list = []

    for i in range(path_num):
        Vt = Vt_list[i]
        file_path = file_path_list[i]

        state_num, mode_num, mode_index_list, same_PES_Nloc, diff_PES_Nloc = read_Nloc_anharmonic_nonadiabatic_coupling(file_path)

        # remove first index, which stands for electronic dof.
        mode_index_list = [ i[1:] for i in mode_index_list ]
        mode_num = mode_num -1

        nonadiabatic_T = []
        for j in range(state_num):
            qn_list = mode_index_list[j]
            nonadiabatic_T_element  = analyze_nonadiabatic_transition_factor(qn_list, Vt)
            nonadiabatic_T.append(nonadiabatic_T_element)

        nonadiabatic_T_list.append(nonadiabatic_T)
        diff_PES_Nloc_list.append(diff_PES_Nloc)

    # transform list to shape [state_num , Vt_num]
    nonadiabatic_T_list = np.transpose(nonadiabatic_T_list)
    diff_PES_Nloc_list = np.transpose(diff_PES_Nloc_list)

    color_list = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

    fig = plt.figure(figsize=(20, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=2, figure=fig)
    ax = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])

    for i in range(5):
        ax.plot(Vt_list, nonadiabatic_T_list[i] , linestyle = 'solid', marker = 'o', linewidth = 2, color = color_list[i] , label = str(mode_index_list[i]))
        ax.plot(Vt_list, diff_PES_Nloc_list[i], linestyle = 'dashed', marker = 'o',  linewidth = 2 ,color = color_list[i] )

    ax.legend(loc = 'best')
    ax.set_xlabel('Vt')
    ax.set_ylabel(" T' or $N_{loc}$' " )

    # plot log(N_{loc}') against log(T') in log scale
    nonadiabatic_T_list_flatten = nonadiabatic_T_list.flatten()
    nonadiabatic_T_list_flatten_log = np.log(nonadiabatic_T_list_flatten)
    diff_PES_Nloc_list_flatten = diff_PES_Nloc_list.flatten()
    diff_PES_Nloc_list_flatten_log = np.log(diff_PES_Nloc_list_flatten)

    ax2.scatter( nonadiabatic_T_list_flatten_log  ,diff_PES_Nloc_list_flatten_log , edgecolors = 'black', facecolors = 'none', s=50 )

    ax2.set_xlabel("log(T')")
    ax2.set_ylabel("log($N_{loc}'$)")

    # pick data points log(T') > -2 and linear fit function about log(N_{loc}) and log(T')
    length = len(nonadiabatic_T_list_flatten)
    log_T_cutoff = -2.5
    index_list = [index for index in range(length) if nonadiabatic_T_list_flatten_log[index] > log_T_cutoff]
    nonadiabatic_T_list_flatten_log_slice = nonadiabatic_T_list_flatten_log [index_list]
    diff_PES_Nloc_list_flatten_log_slice = diff_PES_Nloc_list_flatten_log [index_list]

    z = np.polyfit(nonadiabatic_T_list_flatten_log_slice, diff_PES_Nloc_list_flatten_log_slice, 1)
    a1, a0 = z
    x_range = np.linspace(np.min(nonadiabatic_T_list_flatten_log ) , np.max(nonadiabatic_T_list_flatten_log ) , 100)
    y_range = a0 + a1 * x_range

    ax2.plot(x_range, y_range, linewidth = 2 )

    plt.show()
    if (save_bool):
        file_name = "T(solid) vs Nloc(dashed).svg"
        file_name = os.path.join(folder_path, file_name)
        fig.savefig(file_name)