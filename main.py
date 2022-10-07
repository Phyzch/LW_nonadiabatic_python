from util import *

from Generate_input_file import Generate_input_files
from Analyze_IPR import Analyze_IPR_all_state,Read_IPR_all_state, Read_local_density_of_state , Search_state, compare_IPR_subroutine
from Analyze_Multiple_result import Analyze_multiple_simulation_result
from compare_result_with_different_condition import Analyze_IPR_different_condition
from Overlap_of_displaced_state import plot_Franck_condon_factor , analyze_single_Franck_condon_factor, plot_Bessel_function_envolope,  analyze_coupling_strength
from Overlap_of_displaced_state import compute_LW_factor
from Analyze_survival_prob import Analyze_survival_prob_all_state, compare_survival_prob_subroutine
from Analyze_Nloc_with_dilution_factor import analyze_dilution_factor_and_Nloc, analyze_dilution_factor_and_Nloc_change_V0

def main():
    matplotlib.rcParams.update({'font.size': 20})

    # plot_IPR_and_survival_prob_simulation_result()

    # Analyze_multiple_simulation_result()

    # Analyze_IPR_different_condition()

    # compare_survival_prob_and_IPR()

    # analyze dilution factor with Nloc : for LW transition
    analyze_dilution_factor_and_Nloc()
    # analyze_dilution_factor_and_Nloc_change_V0()

    # input_folder_path ="/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/Duschrinsky_rotation/Batch_simulation_45_long_time/"
    # Generate_input_files(input_folder_path)

    # compute Franck condon factor
    # plot_Franck_condon_factor()
    # analyze_single_Franck_condon_factor()
    # plot_Bessel_function_envolope()
    # analyze_coupling_strength()

    # compute LW factor
    # compute_LW_factor()



    plt.show()


def plot_IPR_and_survival_prob_simulation_result():

    parent_file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/"
    # folder_path = "Bchl 5mode/batch_simulation_dE=600/Vt=500,V0=300,a=0.3/"
    folder_path = "BChl try/"

    folder_path = os.path.join(parent_file_path, folder_path)

    save_bool = False

    # Analyze IPR
    Analyze_IPR_all_state(folder_path, save_bool)
    # Analyze survival prob.
    Analyze_survival_prob_all_state(folder_path, save_bool)

def compare_survival_prob_and_IPR():
    save_bool = True
    folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/Bchl 5mode/Vt=353/"

    # file_path1 = "Vt=0,V0=300,a=0.3"
    file_path1 = "Vt=353,V0=0"
    file_path2 = "Vt=353,V0=300,a=0.2"
    file_path3 = "Vt=353,V0=300,a=0.3"

    file_path_list = [ file_path1, file_path2, file_path3]
    file_path_list = [os.path.join(folder_path, path) for path in file_path_list]


    label1 = 'without IVR coupling'
    label2 =  'with IVR coupling, V=300, a=0.2'
    label3 =  'with IVR coupling, V=300, a=0.3'

    label_list = [ label1, label2, label3]

    mode_index= 2

    compare_survival_prob_subroutine(file_path_list, label_list, mode_index, folder_path, save_bool )
    # compare_IPR_subroutine(file_path_list, label_list, mode_index, folder_path, save_bool)

if __name__ == '__main__':
    main()

