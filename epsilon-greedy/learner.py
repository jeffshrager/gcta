from os import sys
import numpy as np
import random, datetime, time, json
from itertools import combinations as combos, product as iterproduct
from math import floor, ceil, sqrt
from sklearn import linear_model
from sklearn import naive_bayes
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import mean_squared_error

class Learner:
    def __init__(self, settings):
        self.settings = settings
        self.num_markers = len(self.settings.markers)
        self.num_drugs = len(self.settings.drugs)

    def learn(self, data):
        xarrays, yarrays = self.parse_data(data)

        lrs = [linear_model.LinearRegression() for i in range(self.num_drugs)]

        # train the linear regression models
        for i in range(len(lrs)):
            if len(xarrays[i]) > 10:
                lrs[i].fit(xarrays[i], yarrays[i])

        accuracies = []
        for i in range(len(lrs)):
            accuracies.append(self.score_accuracy(lrs[i], xarrays[i], yarrays[i]))
        print("Accuracy:", sum(accuracies)/len(accuracies))
        print(accuracies)

        pseudo_weight_matrix = [[0 for i in range(self.num_drugs)] for j in range(self.num_markers)]
        coeffs = []
        for i in range(len(lrs)):
            if len(xarrays[i]) > 10:
                coeffs.append(lrs[i].coef_)
            else:
                coeffs.append([0 for i in range(self.num_markers)])

        coeffs = [list(row) for row in coeffs]

        for i in range(len(coeffs)):
            row = coeffs[i]
            for j in range(len(row)):
                pseudo_weight_matrix[j][i] = row[j]


        print(pseudo_weight_matrix)
        rms_vector = []
        for i in range(len(lrs)):
            preds = lrs[i].predict(xarrays[i])
            true_res = yarrays[i]
            rms_vector.append(mean_squared_error(preds, true_res))

        return pseudo_weight_matrix, rms_vector

    def score_accuracy(self, lr, x, y):
        predicted = cross_val_predict(lr, x, y, cv=10)
        correct = 0
        j = 0
        for onex in x:
            if y[j] > 0 and predicted[j] > 0 or y[j] < 0 and predicted[j] < 0:
                correct += 1
            j += 1
        accuracy = correct/j
        return accuracy

    def parse_data(self, data):
        # data is a list of dicts where each dict goes like this:
        # update_dict = {'px': patient.pxnum, 'month': patient.month, 'tumor_dec': self.settings.rep_tumor(patient.tumor),
        #                 'drug_combo_num': di+1, 'drug_name': patient.treatment[di], 'px_state': patient.state,
        #                 'k_score_prev': k_score_prev, 'k_score_new': k_score_new, 'delta_k': delta_k}
        xs = [[] for i in range(self.num_drugs)]
        ys = [[] for i in range(self.num_drugs)]

        for d in data:
            tumor = self.settings.decode_tumor(d['tumor_dec'])
            # print('==================================================')
            drug_index = self.settings.drugs[d['drug_name']]['index']
            xs[drug_index].append(tumor)
            ys[drug_index].append(d['delta_k'])

        #turn xs and ys into numpy arrays for sklearn later...
        xarrays = [np.asarray(x) for x in xs]
        yarrays = [np.asarray(y) for y in ys]

        return xarrays, yarrays

class Settings:
    def __init__(self, path_to_params):
        self.params = json.load(open(path_to_params, 'r'))
        self.drugs = self.params['drugs']
        self.drugs_lookup_list = self.gen_lookup_list(self.drugs)
        print(self.drugs_lookup_list)
        self.markers = self.params['markers']
        if self.params['num_additional_markers'] > 0:
            self.gen_markers()
        self.markers_lookup_list = self.gen_lookup_list(self.markers)
        print(self.markers_lookup_list)

    def gen_markers(self):
        labels = []
        for i in range(self.params['num_additional_markers']):
            labels.append('X' + str(i))
        for i in range(len(labels)):
            self.markers[labels[i]] = {"index": i+6, "exponent": 0}
        for marker in self.markers:
            self.markers[marker]['exponent'] = (len(self.markers)-1)-self.markers[marker]['index']

    def gen_lookup_list(self, transform):
        lookup = [None for i in range(len(transform))]
        for item in transform:
            lookup[transform[item]['index']] = item
        return lookup

    def rep_tumor(self, markers):
        dec_rep = 0
        for marker in markers:
            dec_rep += 2**self.markers[marker]['exponent']
        return dec_rep

    def decode_tumor(self, encoding):
        format_string = "{0:0"+str(len(self.markers))+"b}"
        tumor = list(format_string.format(encoding))
        for i in range(len(tumor)):
            tumor[i] = int(tumor[i])
        return tumor


