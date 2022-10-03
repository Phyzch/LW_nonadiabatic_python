import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import special

def Overlap_n_m(alpha,m, n):
    # m in quanta in spin down. n ins quanta in spin up
    # exp(- \alpha^2 / 2) <m| exp(-\Delta alpha b^{+}) exp(\Delta alpha b) |n>
    # <m| alpha; n>
    Minimum_quanta = np.min( (m,n) )
    X = 0
    for i in range(Minimum_quanta + 1):
        x = pow(-1, m - i) * pow(alpha, m+n - 2* i) / ( np.math.factorial(int(m-i))  * np.math.factorial(int(n - i))  ) \
        * np.sqrt( float(np.math.factorial(m)) * np.math.factorial(n) ) / np.math.factorial(i)

        X = X + x

    X = X * np.exp(- pow(alpha,2) / 2 )

    return X

def Overlap_n_m_approx(alpha, m ,n):
    if m < n:
        l = m
        m = n
        n = l

    # factor_approx = special.jv( m - n , 2 * (alpha * np.sqrt(n)) ) / np.power(alpha * np.sqrt(n) , m - n) * np.math.factorial(m-n)
    #
    # result1 = np.exp(- np.power(alpha,2)) * np.power( (np.power(alpha,2) * n) , (m-n)/2) / (np.math.factorial(m-n)) * factor_approx

    result1 = np.exp(-  np.power(alpha,2)) * special.jv(m-n ,2 * (alpha * np.sqrt(n)) )
    return result1

def Overlap_n_m_component(alpha, m, n):
    # m in quanta in spin down. n ins quanta in spin up
    # exp(- \alpha^2 / 2) <m| exp(-\Delta alpha b^{+}) exp(\Delta alpha b) |n>
    # <m| alpha; n>
    Minimum_quanta = np.min((m, n))
    X = 0
    index_list = np.array(range(Minimum_quanta+1)) + 1
    component_list = np.zeros([Minimum_quanta+1])
    for i in range(Minimum_quanta + 1):
        x = pow(-1, m - i) * pow(alpha, m + n - 2 * i) / (np.math.factorial(int(m - i)) * np.math.factorial(int(n - i))) \
            * np.sqrt(float(np.math.factorial(m)) * np.math.factorial(n)) / np.math.factorial(i)

        X = X + x
        component_list[i] = x

    X = X * np.exp(- pow(alpha, 2) / 2)

    # plot component
    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    ax.plot(index_list, component_list , marker = 'o')
    ax.set_xlabel('i')

    return X

def estimation_n_m_alpha(alpha, m , n):
    '''
    estimate franck condon factor using last two terms when alpha is small
    :param alpha:
    :param m:
    :param n:
    :return:
    '''
    if m < n:
        l = m
        m = n
        n = l

    result = np.exp(-np.power(alpha,2)/2) * np.sqrt(np.math.factorial(m) / np.math.factorial(n))/ np.math.factorial(m-n) * np.power(alpha, m-n) \

    factor = (1 - np.power(alpha,2) * n / (m-n+1)
                + np.power(alpha,4) * n * (n-1)/(2 * (m-n+1) * (m-n + 2))
                - np.power(alpha,6) * n * (n-1) * (n-2) / ( 6 * (m-n+1) * (m-n+2) * (m-n+3) ) )

    print("factor: " + str(round(factor,2)))
    result = result * factor

    factor_approx = special.jv( m - n , 2 * (alpha * np.sqrt(n)) ) / np.power(alpha * np.sqrt(n) , m - n) * np.math.factorial(m-n)

    print("factor approx: " + str(round(factor_approx , 2)))

    result1 = np.exp(- np.power(alpha,2)) * np.power( (np.power(alpha,2) * n) , (m-n)/2) / (np.math.factorial(m-n)) * factor_approx


    return result , result1
def analyze_single_Franck_condon_factor():
    alpha = 0.5
    n = 10
    m = 10
    X = Overlap_n_m_component(alpha, m, n)
    print("<n|alphpa;m> result: " + str(round(X,2)))

    X_estimate , X_estimate2 = estimation_n_m_alpha(alpha, m, n)
    print("<n|alpha;m> estimate : " + str(round(X_estimate , 2)))
    print("<n|alpha;m> estimate2 : " + str(round(X_estimate2 , 2)))

    plt.show()
def plot_Franck_condon_factor():
    '''
    compute franck condon factor for <m| alpha;n>
    :return:
    '''
    alpha = 0.5

    qn_cutoff = 15

    franck_condon_table = np.zeros([qn_cutoff , qn_cutoff])
    franck_condon_table_approx = np.zeros([qn_cutoff, qn_cutoff ])

    for i in range(qn_cutoff):
        for j in range(qn_cutoff):
            franck_condon_factor = Overlap_n_m(alpha, i, j)
            franck_condon_table[i,j] = franck_condon_factor

            franck_condon_factor_approx = Overlap_n_m_approx(alpha , i, j)
            franck_condon_table_approx [i,j] = franck_condon_factor_approx

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    cp = ax.imshow(franck_condon_table)
    fig.colorbar(cp)

    col_index_list = np.array(range(qn_cutoff))
    row_index_to_plot = 5


    fig1 = plt.figure(figsize=(10, 10))
    spec1 = gridspec.GridSpec(nrows=1, ncols=1, figure=fig1)
    ax1 = fig1.add_subplot(spec1[0,0])

    ax1.plot(col_index_list , np.abs(franck_condon_table_approx[row_index_to_plot]) , label = 'n = ' + str(row_index_to_plot) + " approx" , marker = 'o')

    ax1.plot( col_index_list, np.abs(franck_condon_table[row_index_to_plot]) , label = 'n = ' + str(row_index_to_plot) , marker = 'o')
    ax1.legend(loc = 'best')

    print( np.round( franck_condon_table, 3) )

    plt.show()

def plot_Bessel_function_envolope():
    '''

    :return:
    '''
    N = 5
    x_range = 3
    x_num = 100
    x = np.linspace(0,x_range,x_num)

    bessel_func_value_list = np.zeros([N, x_num])

    for i in range(N):
        bessel_func_value = special.jv(i, x)
        bessel_func_value_list[i] = bessel_func_value

    fig = plt.figure(figsize=(10, 10))
    spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])

    for i in range(N):
        ax.plot(x, bessel_func_value_list[i] , label = 'n=' + str(i))

    ax.legend(loc = 'best')
    plt.show()

