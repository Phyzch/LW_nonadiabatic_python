'''
plot Nloc (both for anharmonic and nonadiabatic coupling) vs dilution factor
'''
from Analyze_Nloc import read_Nloc_anharmonic_nonadiabatic_coupling
from Analyze_survival_prob import compute_dilution_factor
from util import *


def plot_dilution_factor_and_Nloc_single_state_change_Vt(same_PES_Nloc, diff_PES_Nloc, mode_num, dilution_factor, ax):
    '''
    plot result for single state with different nonadiabatic coupling.
    :param same_PES_Nloc:
    :param diff_PES_Nloc:
    :param mode_num_list:
    :param dilution_factor:
    :return:
    '''
    param_num = len(same_PES_Nloc)
    for i in range(param_num):
        if same_PES_Nloc[i] >= 1:
            print("Nloc(same PES) >=1 for state: " + str(mode_num))
            return

    # T~/ (1-T)
    criteria = diff_PES_Nloc / (1-same_PES_Nloc)

    # ax.scatter(criteria, dilution_factor , label = str(mode_num[1:]))
    ax.plot(criteria, dilution_factor , label = str(mode_num[1:]) , marker = 'o')

def plot_dilution_factor_and_Nloc_change_Vt(mode_number_list, dilution_factor_list, same_PES_Nloc_list, diff_PES_Nloc_list, save_bool , fig_path):
    '''
    for the cases we change nonadiabatic coupling Vt
    :param mode_number_list:
    :param dilution_factor_list:
    :param same_PES_Nloc_list:
    :param diff_PES_Nloc_list:
    :return:
    '''
    state_num = len(mode_number_list)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    spec.update(hspace=0.5, wspace=0.3)
    ax = fig.add_subplot(spec[0, 0])

    for i in range(0,5):
        same_PES_Nloc = same_PES_Nloc_list[i]
        diff_PES_Nloc = diff_PES_Nloc_list[i]
        mode_number = mode_number_list[i]
        dilution_factor = dilution_factor_list[i]
        plot_dilution_factor_and_Nloc_single_state_change_Vt(same_PES_Nloc, diff_PES_Nloc, mode_number, dilution_factor, ax)

    ax.set_xlabel("$N'_{loc}$/(1-$N_{loc}$)")
    ax.set_ylabel('dilution factor')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([pow(10,-3),10])
    ax.set_ylim([5 * pow(10,-3) , 1])

    ax.legend(loc = 'best' , prop = {'size': 10})

    if save_bool:
        fig_name = "LW criteria vs dilution factor.png"
        fig_path1 = os.path.join(fig_path,fig_name)
        fig.savefig(fig_path1)

def plot_same_PES_Nloc_vs_dilution_factor(same_PES_Nloc, dilution_factor_list):
    '''

    :return:
    '''
    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    spec.update(hspace=0.5, wspace=0.3)
    ax = fig.add_subplot(spec[0, 0])

    same_PES_Nloc_no_Vt = same_PES_Nloc[:,0]
    dilution_factor = dilution_factor_list[:,0]

    ax.scatter(same_PES_Nloc_no_Vt, dilution_factor)
    ax.set_xlabel('$N_{loc}$')
    ax.set_ylabel('dilution factor')
    ax.set_title('localized state dilution factor')


def analyze_dilution_factor_and_Nloc():
    '''

    :return:
    '''
    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_dE=600/"

    file_path1 = "Vt=0,V0=300,a=0.3"
    file_path2 = "Vt=10,V0=300,a=0.3"
    file_path3 = "Vt=50,V0=300,a=0.3"
    file_path4 = "Vt=100,V0=300,a=0.3"
    file_path5 = "Vt=200,V0=300,a=0.3"
    file_path6 = "Vt=300,V0=300,a=0.3"
    file_path7 = "Vt=500,V0=300,a=0.3"
    file_path8 = "Vt=1000,V0=300,a=0.3"

    file_path_list = [ file_path1, file_path2, file_path3, file_path4, file_path5, file_path6, file_path7]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    path_num = len(file_path_list)

    mode_number_list = []
    dilution_factor_list = []
    same_PES_Nloc_list = []
    diff_PES_Nloc_list = []

    for i in range(path_num):
        file_path = file_path_list[i]
        state_energy, mode_number_list, dilution_factor = compute_dilution_factor(file_path)
        state_num, mode_num, mode_number_list, same_PES_Nloc, diff_PES_Nloc = read_Nloc_anharmonic_nonadiabatic_coupling(file_path)

        dilution_factor_list.append(dilution_factor)
        same_PES_Nloc_list.append(same_PES_Nloc)
        diff_PES_Nloc_list.append(diff_PES_Nloc)

    # transpose to [state_num , parameter set num]
    dilution_factor_list = np.transpose(dilution_factor_list)
    same_PES_Nloc_list = np.transpose(same_PES_Nloc_list)
    diff_PES_Nloc_list = np.transpose(diff_PES_Nloc_list)

    plot_dilution_factor_and_Nloc_change_Vt(mode_number_list, dilution_factor_list, same_PES_Nloc_list, diff_PES_Nloc_list, save_bool , fig_path = folder_path)

    # plot dilution factor vs N_loc_same_PES when Vt=0
    plot_same_PES_Nloc_vs_dilution_factor(same_PES_Nloc_list, dilution_factor_list)

def plot_dilution_factor_and_Nloc_single_state_change_V0(same_PES_Nloc, diff_PES_Nloc, mode_num, dilution_factor, ax):
    '''
    result for case changing anharmonic coupling in molecules
    plot result for single state with different nonadiabatic coupling.
    :param same_PES_Nloc:
    :param diff_PES_Nloc:
    :param mode_num_list:
    :param dilution_factor:
    :return:
    '''
    # (1-T)/T~
    criteria = (1 - same_PES_Nloc) / diff_PES_Nloc

    ax.scatter(criteria, dilution_factor , label = str(mode_num))
    # ax.plot(criteria, dilution_factor , label = str(mode_num) , marker = 'o')

def plot_dilution_factor_and_Nloc_change_V0(mode_number_list, dilution_factor_list, same_PES_Nloc_list, diff_PES_Nloc_list):
    '''
    for the cases we change nonadiabatic coupling Vt
    :param mode_number_list:
    :param dilution_factor_list:
    :param same_PES_Nloc_list:
    :param diff_PES_Nloc_list:
    :return:
    '''
    state_num = len(mode_number_list)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    spec.update(hspace=0.5, wspace=0.3)
    ax = fig.add_subplot(spec[0, 0])

    for i in range(state_num):
        same_PES_Nloc = same_PES_Nloc_list[i]
        diff_PES_Nloc = diff_PES_Nloc_list[i]
        mode_number = mode_number_list[i]
        dilution_factor = dilution_factor_list[i]
        plot_dilution_factor_and_Nloc_single_state_change_V0(same_PES_Nloc, diff_PES_Nloc, mode_number, dilution_factor, ax)

    ax.set_xlabel("(1-$N_{loc}$)/$N'_{loc}$")
    ax.set_ylabel('dilution factor')
    ax.set_yscale('log')
    ax.set_xlim([0, None])

    ax.legend(loc = 'best' , prop = {'size': 10})


def analyze_dilution_factor_and_Nloc_change_V0():
    '''
    for the case we change nonadiabatic coupling V0
    :return:
    '''
    save_bool = False
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/batch_simulation_dE=600,Vt=300/"
    file_path1 = "V0=20"
    file_path2 = "V0=50"
    file_path3 = "V0=100"
    file_path4 = "V0=200"
    file_path5 = "V0=300"
    file_path6 = "V0=350"
    file_path7 = "V0=400"
    file_path8 = "V0=500"

    file_path_list = [file_path1, file_path2, file_path3, file_path4, file_path5, file_path6, file_path7, file_path8]

    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]

    path_num = len(file_path_list)

    mode_number_list = []
    dilution_factor_list = []
    same_PES_Nloc_list = []
    diff_PES_Nloc_list = []

    for i in range(path_num):
        file_path = file_path_list[i]
        state_energy, mode_number_list, dilution_factor = compute_dilution_factor(file_path)
        state_num, mode_num, mode_number_list, same_PES_Nloc, diff_PES_Nloc = read_Nloc_anharmonic_nonadiabatic_coupling(
            file_path)

        dilution_factor_list.append(dilution_factor)
        same_PES_Nloc_list.append(same_PES_Nloc)
        diff_PES_Nloc_list.append(diff_PES_Nloc)

    # transpose to [state_num , parameter set num]
    dilution_factor_list = np.transpose(dilution_factor_list)
    same_PES_Nloc_list = np.transpose(same_PES_Nloc_list)
    diff_PES_Nloc_list = np.transpose(diff_PES_Nloc_list)

    plot_dilution_factor_and_Nloc_change_V0(mode_number_list, dilution_factor_list, same_PES_Nloc_list, diff_PES_Nloc_list)