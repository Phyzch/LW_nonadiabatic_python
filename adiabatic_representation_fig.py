from util import *
import scipy.linalg as LA

def plot_diabatic_and_adiabatic_PES():
    '''

    :return:
    '''
    save_bool = True
    matplotlib.rcParams.update({'font.size': 20})

    frequency = 890
    alpha = 0.169
    c_coupling = np.sqrt(2 / frequency) * pow(frequency,2) * alpha
    Vt = 364

    # x is in unit of (c_i^{2} / (m omega_i^{2})) : distance between minimal
    xmin = -5
    xmax = 5
    Nx = 100
    x_coordinate = np.linspace(xmin, xmax, Nx)

    diabatic_1 = pow(c_coupling,2)/(2 * pow(frequency,2)) * pow(x_coordinate - 0.5,2)
    diabatic_2 = pow(c_coupling,2)/(2 * pow(frequency,2)) * pow(x_coordinate + 0.5,2)
    lower_adiabatic = []
    higher_adiabatic = []
    for index in range(Nx):
        H = np.array([ [diabatic_1[index], Vt] , [Vt, diabatic_2[index]] ])
        w, v = LA.eig(H)
        w = np.real(w)
        lower_adiabatic_ele = min(w)
        higher_adiabatic_ele = max(w)
        lower_adiabatic.append(lower_adiabatic_ele)
        higher_adiabatic.append(higher_adiabatic_ele)

    lower_adiabatic = np.array(lower_adiabatic)
    higher_adiabatic = np.array(higher_adiabatic)

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    ax.plot(x_coordinate, diabatic_1, linewidth = 5, color = 'blue', linestyle = 'dotted')
    ax.plot(x_coordinate, diabatic_2, linewidth= 5, color='red' , linestyle = 'dotted')

    ax.plot(x_coordinate, lower_adiabatic , linewidth = 2, color = 'black', linestyle = 'solid' )
    ax.plot(x_coordinate, higher_adiabatic, linewidth=2, color='black', linestyle='solid')

    if save_bool :
        folder_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/result 2022.10.06/paper fig/potential energy surface/"
        file_name = "diabatic and adiabatic representation.svg"
        file_name = os.path.join(folder_path, file_name)
        fig.savefig(file_name)

    plt.show()

# plot_diabatic_and_adiabatic_PES()