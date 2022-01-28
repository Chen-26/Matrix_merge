import pandas as pd
import numpy as np
import csv

# 1. Use the relative path instead of the absolute path
# 2. Naming the directory as "_" to connect, don't use space ' '
matrix_s = pd.DataFrame(pd.read_csv("E:\Data analysis\install\S.csv",
                                    header=0, index_col=0))
matrix_c = pd.DataFrame(pd.read_csv("E:\Data analysis\install/C.csv",
                                    header=0, index_col=0))
lambda_1 = 0.6
# def filter_matrix_c_lambda_1(num):
matrix_c_filter_lambda1 = matrix_c.copy()
matrix_c_filter_lambda1[(matrix_c.abs() >= lambda_1)] = 1

s_gene_set, s_patient_set = set(matrix_s.columns), set(matrix_s.index)
c_gene_set, c_patient_set = set(matrix_c_filter_lambda1.columns), set(matrix_c_filter_lambda1.index)

# get the union for row names and col names
gene_union_set = (s_gene_set | c_gene_set)
patient_union_set = (s_patient_set | c_patient_set)
gene_intersect_set = (s_gene_set & c_gene_set)
patient_intersect_set = (s_patient_set & c_patient_set)

union_df = pd.DataFrame((["{0}#{1}".format(row, col) for col in gene_union_set] for row in patient_union_set),
                        index=patient_union_set, columns=gene_union_set)



def apply_track(x):
    patient, gene = tuple(x.split("#"))
    # process the intersect case
    if patient in patient_intersect_set and gene in gene_intersect_set:
        return matrix_s.loc[patient, gene] | matrix_c_filter_lambda1.loc[patient, gene]
    # process the case which the cases are in the s set
    if patient in s_patient_set and gene in s_gene_set:
        return matrix_s.loc[patient, gene]
    # process the case which the cases are in the c set
    if patient in c_patient_set and gene in c_gene_set:
        return matrix_c_filter_lambda1.loc[patient, gene]
    # If the case is neither in S set nor in C set, then just fill in 0
    else:
        return 0


union_df = union_df.applymap(apply_track)


exit(0)

## I needs to rewrite your codes and the following is your code
## This operation can be implemented using pandas in an elegant way