from util import *

def plot_spectral_density_BChl():
    '''
    data from  J. Chem. Phys. 134, 024506 (2011) Table II.
    :return:
    '''
    frequency_list = np.array( [
        84, 167, 183, 191, 214, 239, 256, 345, 368, 388, 407, 423, 442, 473, 506, 565, 587, 623, 684,
        696, 710, 727, 776, 803, 845, 858, 890, 915, 967, 980, 1001, 1019, 1066, 1089, 1105, 1117, 1137,
        1158, 1180, 1190, 1211, 1229, 1252, 1289, 1378, 1466, 1519, 1539, 1648, 1680
    ]) # in unit of cm^{-1}

    reorganization_energy = np.array([
        1.3, 1.4, 1.3, 3.7, 1.0, 1.9, 1.4, 5.6, 2.2, 1.6, 2.1, 1.2, 1.0, 0.8, 0.8, 4.6, 2.3, 3.6, 1.6, 1.8,
        1.3, 19.3, 7.4, 3.4, 2.1, 2.2, 25.3, 4.4, 2.6, 3.0, 4.1, 9.9, 2.7, 2.3, 2.4, 11.5, 4.8, 11.9, 3.0, 3.7,
        3.1, 2.6, 3.2, 11.2, 11.9, 3.1, 2.6, 3.6, 4.2, 4.6
    ]) # \lambda in the table, in unit of cm^{-1}

    spectral_density = frequency_list * reorganization_energy
    Huang_Rhys_factor = reorganization_energy / frequency_list


    # selected vibrational state
    # selected_frequency = np.array( [
    #     890, 727, 345, 1117, 1158
    # ] )
    #
    # selected_reorganization_energy = np.array([
    #     25.3, 19.3, 5.6, 11.5, 11.9
    # ])

    # 12 modes most significantly coupled to electronic dof. See Table III of J. Chem. Phys. 134, 024506 (2011)
    selected_frequency = np.array(
        [191, 345, 565, 727, 776, 890, 1019, 1117, 1137, 1158, 1289, 1378]
    )
    selected_reorganization_energy = np.array(
        [3.7, 5.6, 4.6, 19.3, 7.4, 25.3, 9.9, 11.5, 4.8, 11.9, 11.2, 11.9]
    )

    selected_spectral_density = selected_frequency * selected_reorganization_energy
    selected_Huang_Rhys_factor = selected_reorganization_energy / selected_frequency

    fig = plt.figure(figsize=(20, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=2, figure=fig)
    spec.update(hspace=0.5, wspace=0.3)

    ax = fig.add_subplot(spec[0,0])
    ax.stem( frequency_list, spectral_density ) # for spectral density
    ax.stem( selected_frequency, selected_spectral_density , linefmt = 'C1-' )

    ax.set_xlabel('$\omega_{i}$')
    ax.set_ylabel('J($\omega_{i}$)')

    ax1 = fig.add_subplot(spec[0,1])
    ax1.stem( frequency_list, Huang_Rhys_factor )
    ax1.stem( selected_frequency, selected_Huang_Rhys_factor ,linefmt = 'C1-' )

    ax1.set_xlabel('$\omega_{i}$')
    ax1.set_ylabel('$S_{i}$')

    plt.show()
