import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

from Generate_input_file import Generate_input_files
from Analyze_IPR import Analyze_IPR_all_state,Read_IPR_all_state, Read_local_density_of_state , Search_state
from Analyze_Multiple_result import Analyze_multiple_simulation_result
from compare_result_with_different_condition import Analyze_IPR_different_condition
from Overlap_of_displaced_state import plot_Franck_condon_factor , analyze_single_Franck_condon_factor, plot_Bessel_function_envolope

def main():
    matplotlib.rcParams.update({'font.size': 20})

    plot_simulation_result()

    # Analyze_multiple_simulation_result()

    # Analyze_IPR_different_condition()

    # input_folder_path ="/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/Duschrinsky_rotation/Batch_simulation_45_long_time/"
    # Generate_input_files(input_folder_path)

    # plot_Franck_condon_factor()
    # analyze_single_Franck_condon_factor()
    # plot_Bessel_function_envolope()

def plot_simulation_result():

    parent_file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/"
    folder_path = "PBI1 4mode/only_nonadiabatic_coupling_t=0.5/"

    folder_path = os.path.join(parent_file_path, folder_path)

    Analyze_IPR_all_state(folder_path)



if __name__ == '__main__':
    main()

