from Analyze_dimer_code.Analyze_dimer_survival_prob import plot_dimer_survival_prob_single_file
from Analyze_dimer_code.Analyze_dimer_batch_simulation_survival_prob import plot_dimer_survival_prob_selected_states_batch_simulation , \
    analyze_dilution_factor_and_T_T_prime_phase_diagram , analyze_dilution_factor_with_anharmonic_T_dimer
from Analyze_dimer_code.generate_sampling_state_input_file import generate_sampling_state_info_file

def Analyze_dimer_code_main():
    '''

    :return:
    '''
    # plot_dimer_survival_prob_single_file()

    # plot_dimer_survival_prob_two_file()

    # plot_dimer_survival_prob_selected_states_batch_simulation()

    # plot T, T' phase diagram for dimer
    analyze_dilution_factor_and_T_T_prime_phase_diagram()

    # analyze T vs dilution factor on single PES
    # analyze_dilution_factor_with_anharmonic_T_dimer()

    # generate sampling state for simulation.
    # generate_sampling_state_info_file()
