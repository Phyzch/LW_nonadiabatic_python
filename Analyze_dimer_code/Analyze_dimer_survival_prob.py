from util import *


def Read_dimer_survival_prob(file_path):
    file_name = os.path.join(file_path,"survival_prob.txt")

    mode_number = []
    time_list = []
    survival_prob_list = []
    with open (file_name) as f:
        data = f.read()
        data = re.split('\n', data)
        data = [i for i in data if i != '']
        datalen = len(data)

        line = data[0]
        line = re.split(' ', line)
        state_energy = float(line[0])

        line = data[1]
        line = re.split(' ' , line)
        mode_dimer1 = [float(i) for i in line if i!= '']

        mode_number.append(mode_dimer1)

        line = data[2]
        line = re.split(' ' , line)
        mode_dimer2 = [float(i) for i in line if i!= '']

        mode_number.append(mode_dimer2)

        line_index = 3

        step_num = int ((datalen - line_index)/2 )

        for i in range(step_num):
            line = data[line_index]
            line = re.split(' ', line)
            time1 = float(line[0])
            time_list.append(time1)
            line_index = line_index + 1

            line = data[line_index]
            line = re.split(' ', line)
            survival_prob = float(line[0])
            survival_prob_list.append(survival_prob)
            line_index = line_index + 1

        mode_number = np.array(mode_number)
        survival_prob_list = np.array(survival_prob_list )
        time_list = np.array(time_list)

        return mode_number, state_energy, time_list, survival_prob_list

def plot_dimer_survival_prob_two_file():
    '''

    :return:
    '''
    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/dimer case/Vt=0/"
    mode_number, state_energy, time_list, survival_prob_list = Read_dimer_survival_prob(file_path)
    label1 = "single PES "

    file_path2 = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/dimer case/Vt=353/"
    mode_number2, state_energy2, time_list2, survival_prob_list2 = Read_dimer_survival_prob(file_path2)
    label2 = "with Vt=353 "

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    label = "|v1> = " + str(mode_number[0]) + " |v2> = " + str(mode_number[1])
    ax.plot(time_list, survival_prob_list, label = label1, linewidth = 2)
    ax.plot(time_list2, survival_prob_list2, label= label2, linewidth = 2 )

    ax.set_yscale('log')
    ax.set_ylim([pow(10,-4), 1.1])
    plt.show()

def plot_dimer_survival_prob_single_file():
    '''

    :return:
    '''
    file_path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/BChl try/"
    mode_number, state_energy, time_list, survival_prob_list = Read_dimer_survival_prob(file_path)
    label1 = "single PES "


    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0, 0])

    label = "|v1> = " + str(mode_number[0]) + " |v2> = " + str(mode_number[1])
    ax.plot(time_list, survival_prob_list, label = label1, linewidth = 2)

    ax.set_yscale('log')
    ax.set_ylim([pow(10,-4), 1.1])
    plt.show()


