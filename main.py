# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

from Generate_input_file import Generate_input_files
from Analyze_IPR import Analyze_IPR_all_state,Read_IPR_all_state, Read_local_density_of_state , Search_state
from Analyze_Multiple_result import Analyze_multiple_simulation_result
from compare_result_with_different_condition import Analyze_IPR_different_condition

def main():
    # plot_simulation_result()

    Analyze_multiple_simulation_result()

    # Analyze_IPR_different_condition()

    # input_folder_path ="/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/Duschrinsky_rotation/Batch_simulation_45/"
    # Generate_input_files(input_folder_path)

def plot_simulation_result():

    parent_file_path = "/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/"
    folder_path = "Duschrinsky_rotation/Batch_simulation_45/200/0.18/"
    # folder_path = "2_coordinate/Tunneling/Batch_simulation/10/0.25/"
    # folder_path = "1_coordinate/Tunneling/Batch_simulation_long_time/10/0.25/"

    folder_path = os.path.join(parent_file_path, folder_path)

    Analyze_IPR_all_state(folder_path)



if __name__ == '__main__':
    main()

