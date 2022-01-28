'''import numpy as np
import math
import A

K = input()
Matrix_M = [[A.row_A[A.col_A[0].index(p)][k] for k in range(0, K + 1)] for p in A.col_A[0]]

a = list()
N = []
V = []
for i in A.col_A[0][range(1, len(A.col_A[0] + 1))]:
    for p in A:
        a.append(A.row_A[[A.row_A[A.col_A[0].index(i)]]][[A.row_A[A.col_A[0].index(p)]]])
        N.append(max(a))
        if N[i] < 0.5:
            b = np.std(A.row_A[i + 1][1:]) / 2 * (math.sqrt(K)) * (np.mean(A.row_A[i][1:]))
            V.append(b)
        else:
            b = np.std(A.row_A[i + 1][1:]) / (np.mean(A.row_A[i][1:]))
            V.append(b)

Co = sum(N)
Me = sum(V)
Weight = Co + Me
print(Weight)
'''

# above is my original code
# i will rewrite them
import random
import copy
import numpy as np

import pandas as pd

matrix_A = pd.DataFrame(pd.read_csv("E:\Data analysis\install/C.csv",
                                    header=0, index_col=0))
# K = input()
K = 5
'''
# a submatrix of matrix_A which is determined by K
m_gene_set = random.sample(matrix_A.columns, K)
m_patient_set = set(matrix_A.index)


matrix_M = pd.DataFrame((["{0}#{1}".format(row, col) for col in m_gene_set] for row in m_patient_set),
                        index=m_patient_set, columns=m_gene_set)
'''

# A more convenient way
sqrt = np.sqrt(K)
matrix_M = matrix_A[, 0: K]
Eta = []
Upsilon = []
Sigma = []
Mu = []
i = 0
# index属性遍历
for index in matrix_M.index:
    maximum = max(matrix_M.loc(index))
    standard = matrix_M.loc(index).std()
    average = matrix_M.loc(index).mean()
    Eta.append(maximum)
    Sigma.append(standard)
    Mu.append(average)

for i in len(Eta):
    if Eta[i] < 0.5:
        value = Sigma[i]/(2 * sqrt * Mu[i])
        Upsilon.append(value)
    else:
        value = Sigma(i) / Mu(i)
        Upsilon.append(value)

Coverage = sum(Eta)
Mutual_exclusivity = sum(Upsilon)
Weight = Coverage + Mutual_exclusivity