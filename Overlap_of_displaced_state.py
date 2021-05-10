import numpy as np

def Overlap_n_m(alpha,m, n):
    # m in quanta in spin down. n ins quanta in spin up
    # exp(- \alpha^2 / 2) <m| exp(-\Delta alpha b^{+}) exp(\Delta alpha b) |n>
    Minimum_quanta = np.min( (m,n) )
    X = 0
    for i in range(Minimum_quanta + 1):
        x = pow(-1, m - i) * pow(alpha, m+n - 2* i) / ( np.math.factorial(int(m-i))  * np.math.factorial(int(n - i))  ) \
        * np.sqrt( np.math.factorial(m) * np.math.factorial(n) ) / np.math.factorial(i)

        X = X + x

    X = X * np.exp(- pow(alpha,2) / 2 )

    return X


# lambda_up = 2000
# lambda_down = 2000
# frequency = 1149
# alpha = (lambda_up + lambda_down) / frequency
# overlap = Overlap_n_m(alpha, 3 ,3)
# print(overlap)
# m = 9
# X_list = []
# n_list = range(0,15)
# for n in n_list:
#     X = Overlap_n_m(alpha, m, n)
#     X_list.append(X)
#
#
# print('alpha ' + str(alpha))
# print(n_list)
# print(X_list)