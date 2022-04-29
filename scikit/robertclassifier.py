import os
import sys
from sklearn.model_selection import cross_val_predict
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
import numpy as np

def main(path_to_sim_output):
    f = open(path_to_sim_output, 'r')
    lines = f.readlines()
    if len(lines) > 30000:
        lines = lines[-30000:]

    Xs = [[] for i in range(11)]
    Ys = [[] for i in range(11)]

    drugs = ["BEVACIZUMAB", "TEMOZOLOMIDE", "CABAZITAXEL", "TCAR", "DISATINIB" , "NIVOLUMAB" , "DOXORUBICIN", "DURVALUMAB", "PEMBROLIZUMAB", "VARLILUMAB",  "SORAFENIB"]

    for line in lines:
        words = line.split()
        tumor_id = int(words[2])
        num = words[3]
        drug = words[4]
        delta_ks = int(words[8])

        #format the tumor type as a binary number/array for input
        tumor_list = list("{0:06b}".format(tumor_id))
        #convert list of chars to list of ints
        for i in range(len(tumor_list)):
            tumor_list[i] = int(tumor_list[i])

        drug_index = drugs.index(drug)

        Xs[drug_index].append(tumor_list)
        Ys[drug_index].append(delta_ks)

    xarrays = [np.asarray(x) for x in Xs]
    yarrays = [np.asarray(y) for y in Ys]

    lrs = [linear_model.LinearRegression() for i in range(11)]

    for i in range(11):
        lrs[i].fit(xarrays[i], yarrays[i])

    pseudo_weight_matrix = [[0 for i in range(11)] for j in range(6)]
    coeffs = [lr.coef_ for lr in lrs]
    coeffs = [list(row) for row in coeffs]

    for i in range(len(coeffs)):
        row = coeffs[i]
        for j in range(len(row)):
            pseudo_weight_matrix[5-j][i] = row[j]

    return pseudo_weight_matrix

if __name__ == '__main__':
    filepath = sys.argv[-1]
    print(main(filepath))
