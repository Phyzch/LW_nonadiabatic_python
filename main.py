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

def main():
    # plot_simulation_result()

    Analyze_multiple_simulation_result()

    # input_folder_path ="/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/Batch_simulation/"
    # Generate_input_files(input_folder_path)

def plot_simulation_result():

    parent_file_path = "/home/phyzch/CLionProjects/4_point_correlation_calculation/result/spin_boson_LW_model/"
    folder_path = "resonant(new)/scaling_0.2_no_electronic_state_coupling/"

    folder_path = os.path.join(parent_file_path, folder_path)

    Analyze_IPR_all_state(folder_path)



if __name__ == '__main__':
    main()

