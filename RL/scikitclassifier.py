from os import sys
import numpy as np
import random
import datetime
import time
import json
from itertools import combinations as combos, product as iterproduct
from math import floor, ceil, sqrt
from sklearn.model_selection import cross_val_predict
from sklearn import linear_model
from sklearn.metrics import mean_squared_error

class Learner:
    def __init__(self, settings, data_file):
        self.settings = settings
        self.data_file = data_file
        self.num_markers = len(self.settings.markers_dict)
        self.num_drugs = len(self.settings.drug_dict)

    def generate_hypotheses(self):
        #get the best matrix
        weight_matrix = self.settings.param('best-matrix')
        reversed_drug_dict = dict((v,k) for k,v in self.settings.drug_dict.items())
        mut_dict = self.settings.markers_dict

        muts = []
        for i in range(len(mut_dict)):
            mut = [0 for j in range(len(mut_dict))]
            mut[i] = 1
            muts.append(mut)

        hypotheses = {}
        for i in range(len(muts)):
            mut = muts[i]
            scores = list(np.matmul(mut, weight_matrix))

            drug_scores = {}
            for j in range(len(scores)):
                drug_scores[reversed_drug_dict[j]] = scores[j]

            hypotheses[mut_dict[i]] = {'drug-scores': drug_scores, 'mean': np.mean(scores), 'stddev': np.std(scores)}
            
        return hypotheses

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

    def run(self, settings_obj_num):
        self.settings.params['ts'] = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

        xarrays, yarrays = self.parse_log()

        meta_results = []
        for k in range(self.settings.param("trial-n")):

            lrs = [linear_model.LinearRegression() for i in range(self.num_drugs)]

            #train the linear regression models
            for i in range(len(lrs)):
                if len(xarrays[i]) > 10:
                    lrs[i].fit(xarrays[i], yarrays[i])

            accuracies = []
            for i in range(len(lrs)):
                if len(xarrays[i]) > 10:
                    x = xarrays[i]
                    y = yarrays[i]
                    lr = lrs[i]
                    accuracies.append(self.score_accuracy(lr, x, y))
            mean_accuracy = np.mean(accuracies)

            meta_results.append(mean_accuracy)
            print("On trial", str(k+1) + '/' + str(self.settings.param('trial-n')), "of experiment configuration", str(settings_obj_num) + ". Accuracy", str(mean_accuracy))

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

            if mean_accuracy > self.settings.param("best-accuracy"):
                self.settings.params['best-accuracy'] = mean_accuracy
                #combine with the confidence, somehow.... then save it
                best_matrix = self.calculate_matrix(pseudo_weight_matrix)
                self.settings.params['best-matrix'] = best_matrix

        self.settings.params['results'] = meta_results
        self.settings.params['mean-results'] = np.mean(meta_results)
        print('Mean results from this configuration:', str(self.settings.params['mean-results']), ' std =',str(np.std(meta_results)))
        self.settings.hypotheses = self.generate_hypotheses()

        return self.settings.params, self.settings.hypotheses
        # return self.settings.params

    def calculate_delta_k_thresholds(self, split_lines):
        """
        Calculates the delta k thresholds and sets them accordingly
        """  
        # Notch finter: Figure the mean and std of the delta-ks, and then move 1 (or 
        # per param: delta-k-range-for-thresholds) each way, and filter out anything between those.
        all_dk = [int(line[-1]) for line in split_lines]

        mean_delta_k = np.mean(all_dk)
        std_delta_k = np.std(all_dk)
        print('mean-delta-k = ',str(mean_delta_k))
        print('std-delta-k = ',str(std_delta_k))

        adjusted_std_delta_k = std_delta_k * self.settings.param('delta-k-std-multiplier-threshold')
        print('adjusted-std-delta-k = ',str(adjusted_std_delta_k))

        self.settings.params['max-delta-k'] = max(all_dk)
        self.settings.params['min-delta-k'] = min(all_dk)

        self.settings.params['delta-k-positive-threshold'] = mean_delta_k + adjusted_std_delta_k
        self.settings.params['delta-k-negative-threshold'] = mean_delta_k - adjusted_std_delta_k

        print('delta-k-range', str(self.settings.param('min-delta-k')), ',', 
              str(self.settings.param('max-delta-k')), )
        print('delta-k-thresholds:', str(self.settings.param('delta-k-negative-threshold')), ',', 
              str(self.settings.param('delta-k-positive-threshold')))

    def parse_line(self, line):
        """
        Parse an individual line and create a dict as such:

        {'tumor': binary representation of the tumor as a list,
        'drug_index': the index of the drug in self.settings.drug_dict,
        'delta_k': the delta_k from the line}
        """
        tumor_id = int(line[2])
        drug_name = line[4]
        delta_k = int(line[-1])

        parsed_line = {}

        format_string = "{0:0"+str(self.num_markers)+"b}"
        x = list(format_string.format(tumor_id))
        #convert list of chars to list of ints
        for i in range(len(x)):
            x[i] = int(x[i])
        parsed_line['tumor'] = x

        #append the drug index for later sorting...
        parsed_line['drug_index'] = self.settings.drug_dict[drug_name]

        #append the delta_k for later learning
        parsed_line['delta_k'] = delta_k

        return parsed_line

    def parse_log(self):
        self.data_file.seek(0)
        lines = self.data_file.readlines()
        split_lines = [line.split() for line in lines]

        # figure out where the last batch of data starts by looking for the last
        # patient number 1
        pt_num = int(split_lines[0][0])
        last_batch_start_index = 0

        for i in range(1, len(split_lines)):
            if int(split_lines[i][0]) == 1 and pt_num != 1:
                last_batch_start_index = i
            pt_num = int(split_lines[i][0])

        # only take the last batch of data
        split_lines = split_lines[last_batch_start_index:]

        #CALCULATE ALL THE STATISTICS
        self.split_lines = [line[:] for line in split_lines]
        # self.calculate_confidence_stats(split_lines)

        self.calculate_delta_k_thresholds(split_lines)
 
        #list of dicts with keys "tumor" (binary rep), "drug_index", "delta_k"
        data_chunks = [self.parse_line(line) for line in split_lines]
 
        xs = [[] for i in range(self.num_drugs)]
        ys = [[] for i in range(self.num_drugs)]

        for chunk in data_chunks:
            delta_k = chunk['delta_k']
            if (delta_k >= self.settings.param('delta-k-positive-threshold') 
                or delta_k <= self.settings.param('delta-k-negative-threshold')):
                tumor_id = chunk['tumor']
                drug_index = chunk['drug_index']
                xs[drug_index].append(tumor_id)
                ys[drug_index].append(delta_k)

        #turn xs and ys into numpy arrays for sklearn later...
        xarrays = [np.asarray(x) for x in xs]
        yarrays = [np.asarray(y) for y in ys]

        print('==========')
        print('xs and ys')
        xs_lens = [len(x) for x in xs]
        ys_lens = [len(y) for y in ys]
        print(xs_lens)
        print(ys_lens)
        print('==========')

        return xarrays, yarrays

    def check_significant(self, v1, e1, v2, e2):
        return (v1>v2 and v1-e1>v2+e2) or (v1<v2 and v1+e1<v2-e2)

    def calculate_matrix(self, pseudo_weight_matrix):
        drugs = [0 for i in range(len(self.settings.drug_dict))]
        for drug in self.settings.drug_dict:
            ind = self.settings.drug_dict[drug]
            drugs[ind] = drug

        markers = [self.settings.markers_dict[i] for i in self.settings.markers_dict]

        stats_dict = {marker: {drug: {'dks': [], 'learner': 0, 'n': 0, 'avg': 0, 'std': 0, 'stderr': 0} for drug in drugs} for marker in markers}

        for line in self.split_lines:
            tumor = int(line[2])
            drug = line[4]
            delta_k = int(line[-1])

            format_string = "{0:0" + str(len(markers)) + "b}"
            tumor_bin = list(format_string.format(tumor))
            tumor_bin = [int(i) for i in tumor_bin]

            mutations = [markers[i] for i in range(len(markers)) if tumor_bin[i]]
            for mut in mutations:
                stats_dict[mut][drug]['dks'].append(delta_k)

        for m in markers:
            for d in drugs:
                stats_dict[m][d]['n'] = len(stats_dict[m][d]['dks'])
                stats_dict[m][d]['avg'] = np.mean(stats_dict[m][d]['dks'])
                stats_dict[m][d]['std'] = np.std(stats_dict[m][d]['dks'])
                stats_dict[m][d]['stderr'] = stats_dict[m][d]['std']/sqrt(stats_dict[m][d]['n'])

        for i in range(len(pseudo_weight_matrix)): #6 markers
            for j in range(len(pseudo_weight_matrix[i])): #11 drugs
                stats_dict[markers[i]][drugs[j]]['learner'] = pseudo_weight_matrix[i][j]

        for m in markers:
            lrnr_stderr_drug_tuples = [(drug, stats_dict[m][drug]['learner'], stats_dict[m][drug]['stderr']) for drug in drugs]
            lrnr_stderr_drug_tuples.sort(key=lambda x: x[1])
            lrnr_stderr_drug_tuples.reverse()

            best_drug, best_lrnr, best_stderr = lrnr_stderr_drug_tuples[0]
            good_drugs = [best_drug]
            boundary_reached = False
            for i in range(1, len(lrnr_stderr_drug_tuples)):
                if boundary_reached:
                    break
                drug, lrnr, stderr = lrnr_stderr_drug_tuples[i]
                if self.check_significant(best_lrnr, best_stderr, lrnr, stderr):
                    boundary_reached = True
                else:
                    good_drugs.append(drug)

            stats_dict[m]['significance-group'] = good_drugs

            print(m, ':', ', '.join(stats_dict[m]['significance-group']))

        min_weight = pseudo_weight_matrix[0][0]
        max_weight = pseudo_weight_matrix[0][0]
        for row in pseudo_weight_matrix:
            if min(row) < min_weight:
                min_weight = min(row)
            if max(row) > max_weight:
                max_weight = max(row)

        min_weight = 0
        max_weight = 1

        confidence_matrix = [[min_weight for j in range(len(drugs))] for i in range(len(markers))]
        for i in range(len(markers)):
            for j in range(len(drugs)):
                if drugs[j] in stats_dict[markers[i]]['significance-group']:
                    confidence_matrix[i][j] = max_weight

        # average confidence_matrix with best_matrix
        # best_matrix = [[(confidence_matrix[i][j]+pseudo_weight_matrix[i][j])/2 for j in range(len(drugs))] for i in range(len(markers))]
        # print('MATRIX BULLSHIT')
        # print('pseudo')
        # print(pseudo_weight_matrix)
        # print('confidence')
        # print(confidence_matrix)
        # print('averaged...')
        # print(best_matrix)
        # return best_matrix

        return confidence_matrix

class ExpSettings:
    def __init__(self, keys, values, drug_dict, markers_dict):
        self.params = dict(zip(keys, values))
        self.drug_dict = drug_dict
        self.markers_dict = markers_dict
        self.params['best-accuracy'] = -1
        self.params['best-matrix'] = []

    def param(self, key):
        return self.params[key]

class Framework:
    def __init__(self, path_to_param_file, timestamp):
        self.timestamp = timestamp
        self.params = json.load(open(path_to_param_file, 'r'))
        self.learner_params = self.params['learner']
        self.settings_objs = []

        self.run_label = self.params['run-label']
        self.log_name = self.params['classifier-output']
        self.data_name = self.params['sim-output']

        self.path_to_run_log = 'experiments/'+self.timestamp+'/'+self.run_label+'/'+self.log_name
        self.run_log = open(self.path_to_run_log, 'a')
        self.path_to_data_file = 'experiments/'+self.timestamp+'/'+self.run_label+'/'+self.data_name
        self.data_file = open(self.path_to_data_file, 'r')

        self.drug_dict_raw = self.params['simulator']['drugs']
        self.drug_dict = self.make_drug_dict() #key: drug, value: index

        self.markers_raw = self.params['simulator']['markers'] #list of markers
        self.markers_dict = self.make_markers_dict() #key: index, value: marker

    def make_drug_dict(self):
        """
        Let these drugs never be out of order again ffs
        """
        drug_dict = {}
        for drug in self.drug_dict_raw:
            drug_dict[drug.upper()] = self.drug_dict_raw[drug]['index']
        return drug_dict

    def make_markers_dict(self):
        """
        For internal representation of markers to numbers in array...
        """
        markers_dict = {}
        for i in range(len(self.markers_raw)):
            markers_dict[i] = self.markers_raw[i]
        return markers_dict

    def pnl(self, to_pnl):
        self.run_log.write(to_pnl+'\n')
        print(to_pnl)

    def create_settings_objs(self):
        keys = list(self.learner_params.keys())
        raw_vals = list(self.learner_params.values())
        # params_products = iterproduct(*raw_vals)

        #remove listy bits
        self.settings_objs = [ExpSettings(keys, raw_vals, self.drug_dict, self.markers_dict)]

        # for prod in params_products:
        #     vals = [v for v in prod]
        #     self.settings_objs.append(ExpSettings(keys, vals, self.drug_dict, self.markers_dict))

    def run_experiments(self):
        for i in range(len(self.settings_objs)):
            settings_obj = self.settings_objs[i]
            l = Learner(settings_obj, self.data_file)
            to_log, hypotheses = l.run(i+1)
            # to_log = l.run(i+1)
            # hf = open('experiments/'+self.run_label+'/hypotheses-outputs', 'a')
            # hf.write(str(hypotheses).replace('\'', '"').upper()+'\n\n')
            # hf.close()
            print(hypotheses)
            to_log['drug-dict'] = settings_obj.drug_dict
            self.log_run(to_log)
        self.cleanup()
        return 'Classifier done!'

    def log_run(self, to_log):
        self.run_log.write(str(to_log).replace('\'', '"'))
        self.run_log.write('\n')

    def cleanup(self):
        self.run_log.close()

def main(path_to_params, timestamp):
    frame = Framework(path_to_params, timestamp)
    frame.create_settings_objs()
    return frame.run_experiments()

if __name__ == "__main__":
    main('parameters.json')
