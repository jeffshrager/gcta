from os import sys
import tensorflow as tf
import numpy as np
import random
import datetime
import time
import json
from itertools import combinations as combos, product as iterproduct
from math import floor, ceil

class Learner:
    def __init__(self, settings, data_file):
        self.settings = settings
        self.data_file = data_file
        self.num_markers = len(self.settings.markers_dict)
        self.num_drugs = len(self.settings.drug_dict)

    def generate_tf_variable(self, shape, mode):
        if mode == 'zeros':
            return tf.Variable(tf.zeros(shape))
        elif mode == 'noisy':
            return tf.Variable(tf.truncated_normal(shape, stddev=0.1))
        else:
            raise Exception("Bad mode for tf variable generation")

    def generate_hypotheses(self):
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
            scores_copy = scores[:]
            scores_copy.sort()

            drug_scores = {}
            for j in range(len(scores)):
                drug_scores[reversed_drug_dict[j]] = scores[j]

            hypotheses[mut_dict[i]] = {'drug-scores': drug_scores, 'mean': np.mean(scores), 'stddev': np.std(scores)}
            
        return hypotheses

    def run(self, settings_obj_num):
        self.settings.params['ts'] = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

        training_x, training_y, testing_x, testing_y = self.parse_log()
        self.settings.params['training-length'] = len(training_x)
        self.settings.params['testing-length'] = len(testing_x)

        num_training_samples = self.settings.param('num-data-to-learn')
        # num_training_samples = int(round(len_training_x * self.settings.param("fraction-training-samples"))) #removed param
        print('num_training_samples',num_training_samples)

        meta_results = []
        for k in range(self.settings.param("trial-n")):
            g = tf.Graph()

            with g.as_default():
                tf.set_random_seed(random.randint(0, 10000000))

                sess = tf.InteractiveSession()    

                x = tf.placeholder(tf.float32, [None, self.num_markers]) #input (tumors)

                num_hidden_units = self.settings.param("num-hidden-units")

                W1 = self.generate_tf_variable([self.num_markers, num_hidden_units], 'noisy')
                b1 = self.generate_tf_variable([num_hidden_units], 'noisy')

                h = tf.matmul(x, W1) + b1

                W2 = self.generate_tf_variable([num_hidden_units, self.num_drugs], 'noisy')
                b2 = self.generate_tf_variable([self.num_drugs], 'noisy')

                y = tf.matmul(h, W2) + b2

                y_ = tf.placeholder(tf.float32, [None, self.num_drugs])
                cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y))

                train_step = tf.train.GradientDescentOptimizer(self.settings.param("gd-step-size")).minimize(cross_entropy)
                
                # check python version and run according initializer for tf
                if sys.version[0] == '2':
                    tf.initialize_all_variables().run()
                else:
                    tf.global_variables_initializer().run()

                # Train
                for i in range(self.settings.param("num-training-cycles")):
                    #choose num_training_samples random points from training_x and training_y...
                    sampled = random.sample(range(1, len(training_x)), num_training_samples)
                    batch_xs = [training_x[i] for i in sampled]
                    batch_ys = [training_y[i] for i in sampled]
                    sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})

                # Test 
                max_y = tf.nn.top_k(y, 2).indices #[+, 0, 0, -, 0, -, 0, 0, +, 0, 0]
                max_y_ = tf.nn.top_k(y_, 2).indices
                max_max_y = tf.nn.top_k(max_y, 2).values
                max_max_y_ = tf.nn.top_k(max_y_, 2).values
                correct_prediction = tf.equal(max_max_y, max_max_y_)

                # correct_prediction = tf.equal(tf.nn.top_k(y, 2).indices, tf.nn.top_k(y_, 2).indices)
                accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
                result = sess.run(accuracy, feed_dict={x: testing_x, y_: testing_y}) #[0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0]

                # Record
                meta_results.append(result)
                if meta_results[-1] > self.settings.param("best-accuracy"):
                    self.settings.params['best-accuracy'] = meta_results[-1]
                    best_matrix = sess.run(W)
                    self.settings.params['best-matrix'] = [list(ar) for ar in list(best_matrix)]
                print('On trial', str(k+1) + '/' + str(self.settings.param('trial-n')), 'of experiment configuration', str(settings_obj_num) + '. Accuracy:', str(result))

        self.settings.params['results'] = meta_results
        self.settings.params['mean-results'] = np.mean(meta_results)
        print('Mean results from this configuration:', str(self.settings.params['mean-results']), ' std=',str(np.std(meta_results)))
        self.settings.hypotheses = self.generate_hypotheses()

        return self.settings.params, self.settings.hypotheses

    def parse_log(self):
        self.settings.params['max-delta-k'] = 0
        self.settings.params['min-delta-k'] = 0
        training_x = []
        testing_x = []
        training_y = []
        testing_y = []

        self.data_file.seek(0)
        lines = self.data_file.readlines()
        split_lines = [line.split() for line in lines]

        # figure out where the last batch of data starts by looking for the last
        # patient number 1
        pt_num = int(split_lines[0][0])
        last_batch_start_index = 0

        for i in range(1, len(split_lines)):
            if int(split_lines[i][0]) == 1 and pt_num != 1:
                print('FOUND THE BREAK')
                last_batch_start_index = i
            pt_num = int(split_lines[i][0])

        print('last batch index', last_batch_start_index)

        # only take the last batch of data
        split_lines = split_lines[last_batch_start_index:]

        print('len of split_lines', len(split_lines))

        #index 3 is the combo number...
        # data_chunks = [[line] for line in split_lines]
        data_chunks = []
        current_chunk = [split_lines[0]]
        for line in split_lines[1:]:
            if int(line[3]) == 1:
                data_chunks.append(current_chunk[:])
                current_chunk = [line]
            else:
                current_chunk.append(line)

        # Notch finter: Figure the mean and std of the delta-ks, and then move 1 (or 
        # per param: delta-k-range-for-thresholds) each way, and filter out anything between those.

        all_dk = [int(chunk[0][-1]) for chunk in data_chunks]
        mean_delta_k = np.mean(all_dk)
        std_delta_k = np.std(all_dk)
        print('mean-delta-k = ',str(mean_delta_k))
        print('std-delta-k = ',str(std_delta_k))
        adjusted_std_delta_k = std_delta_k * self.settings.param('delta-k-std-multiplier-threshold')
        print('adjusted-std-delta-k = ',str(adjusted_std_delta_k))
        self.settings.params['max-delta-k'] = max(all_dk)
        self.settings.params['min-delta-k'] = min(all_dk)
        # !!!!!!!!!!! Trying an experiment forcing this to zero mean:
        # mean_delta_k = 0
        self.settings.params['delta-k-positive-threshold'] = mean_delta_k + adjusted_std_delta_k
        self.settings.params['delta-k-negative-threshold'] = mean_delta_k - adjusted_std_delta_k

        print('len(data_chunks)',len(data_chunks))
        test_training_break_point = int(len(data_chunks)*(1-self.settings.param('test-fraction')))
        print('test_training_break_point', test_training_break_point)
        training_data = data_chunks[:test_training_break_point]
        testing_data = data_chunks[test_training_break_point:]

        print('delta-k-range', str(self.settings.param('min-delta-k')), ',', 
              str(self.settings.param('max-delta-k')), )
        print('delta-k-thresholds:', str(self.settings.param('delta-k-negative-threshold')), ',', 
              str(self.settings.param('delta-k-positive-threshold')))

        #loop through training chunks to produce training_x, training_y
        for chunk in training_data:
            x, y_, delta_k = self.parse_chunk(chunk)
            if (delta_k >= self.settings.param('delta-k-positive-threshold') 
                or delta_k <= self.settings.param('delta-k-negative-threshold')):
               training_x.append(x)
               training_y.append(y_)

        #loop through testing chunks to produce testing_x, testing_y
        for chunk in testing_data:
            x, y_, delta_k = self.parse_chunk(chunk)
            if delta_k >= self.settings.param('delta-k-positive-threshold'):
                testing_x.append(x)
                testing_y.append(y_)

        print('len of testing chunks:', str(len(testing_x)))

        if (len(testing_x) == 0):
            raise Exception("Too few testing samples to run!")

        return training_x, training_y, testing_x, testing_y

    def parse_chunk(self, chunk):
        #these two values should be the same in all given lines
        tumor_id = int(chunk[0][2])
        delta_k = int(chunk[0][-1])

        drug_names = [line[4] for line in chunk]

        #format the tumor type as a binary number/array for input
        format_string = "{0:0"+str(self.num_markers)+"b}"
        x = list(format_string.format(tumor_id))
        #convert list of chars to list of ints
        for i in range(len(x)):
            x[i] = int(x[i])

        if self.settings.param("background-mode") == 'random':
            y_ = [((random.random()*self.settings.param("random-background-multiplier"))-self.settings.param("random-background-shift")) for i in range(self.num_drugs)]
        elif self.settings.param("background-mode") == 'constant':
            y_ = [self.settings.param("background-constant") for i in range(self.num_drugs)]
        else:
            raise Exception("Bad background-mode :(")

        if self.settings.param("training-mode") == 'constant':
            if delta_k > 0:
                for drug_name in drug_names:
                    y_[self.settings.drug_dict[drug_name]] = self.settings.param('training-constant-positive')
            elif delta_k < 0:
                for drug_name in drug_names:
                    y_[self.settings.drug_dict[drug_name]] = self.settings.param('training-constant-negative')
            else:
                for drug_name in drug_names:
                    y_[self.settings.drug_dict[drug_name]] = 0.0
        elif self.settings.param('training-mode') == 'random':
            for drug_name in drug_names:
                y_[self.settings.drug_dict[drug_name]] = ((random.random()*self.settings.param("random-training-multiplier"))-self.settings.param("random-training-shift"))
        elif (self.settings.param("training-mode") == 'true-random'):
            #overwrites whatever default values are set with self.settings.param("background_mode")
            #replaces with random noise for all values
            y_ = [(random.random()*self.settings.param('random-training-multiplier'))-self.settings.param('random-training-shift') for i in range(self.num_drugs)]
        elif (self.settings.param("training-mode") == 'normal'):
            # "normalized" delta k (i.e. divide by delta_k scale size)
            # this is bullshit
            for drug_name in drug_names:
                y_[self.settings.drug_dict[drug_name]] = delta_k/self.settings.param('normal-training-divisor')
        elif (self.settings.param("training-mode") == 'stretch'):
            for drug in drug_names:
                y_[self.settings.drug_dict[drug_name]] = (delta_k+self.settings.param('stretch-training-shift'))/self.settings.param('stretch-training-divisor')
        else:
            raise Exception("Bad training-mode")

        return x, y_, delta_k

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
    def __init__(self, path_to_param_file):
        self.ts = round(time.time())
        self.params = json.load(open(path_to_param_file, 'r'))
        self.learner_params = self.params['learner']
        self.settings_objs = []
        self.run_label = self.params['run-label']
        self.log_name = self.params['classifier-output']
        self.data_name = self.params['sim-output']
        # self.data_name = 'last-quarter.xls' #manual entry for which data file to use
        self.run_log = open('experiments/'+self.run_label+'/'+self.log_name, 'a')
        self.path_to_data_file = 'experiments/'+self.run_label+'/'+self.data_name
        print('data file: ' + self.path_to_data_file)
        self.data_file = open(self.path_to_data_file)
        self.drugs = list(self.params['simulator']['drugs'].keys())
        self.drugs = [drug.upper() for drug in self.drugs]
        self.drug_dict = dict(zip(self.drugs, range(len(self.drugs))))
        print('drug dict')
        print(self.drug_dict)
        self.markers = self.params['simulator']['markers'] #list of markers
        self.markers_dict = self.make_markers_dict()

    def make_markers_dict(self):
        """
        For internal representation of markers to numbers in array...
        """
        markers_dict = {}
        for i in range(len(self.markers)):
            markers_dict[i] = self.markers[i]
        return markers_dict

    def pnl(self, to_pnl):
        self.run_log.write(to_pnl+'\n')
        print(to_pnl)

    def create_settings_objs(self):
        keys = list(self.learner_params.keys())
        raw_vals = list(self.learner_params.values())
        # params_products = iterproduct(*raw_vals)

        # remove listy bits
        self.settings_objs = [ExpSettings(keys, raw_vals, self.drug_dict, self.markers_dict)]

        # for prod in params_products:
        #     vals = [v for v in prod]
        #     self.settings_objs.append(ExpSettings(keys, vals, self.drug_dict, self.markers_dict))

    def run_experiments(self):
        for i in range(len(self.settings_objs)):
            settings_obj = self.settings_objs[i]
            l = Learner(settings_obj, self.data_file)
            to_log, hypotheses = l.run(i+1)
            hf = open('experiments/'+self.run_label+'/hypotheses-outputs', 'a')
            hf.write(str(hypotheses).replace('\'', '"').upper()+'\n\n')
            hf.close()
            # to_log['hypotheses'] = str(hypotheses).replace('\'', '"').upper()
            to_log['drug-dict'] = settings_obj.drug_dict
            self.log_run(to_log)
        self.cleanup()
        return 'Classifier done!'

    def log_run(self, to_log):
        self.run_log.write(str(to_log).replace('\'', '"'))
        self.run_log.write('\n')

    def cleanup(self):
        self.run_log.close()

def main(path_to_params):
    frame = Framework(path_to_params)
    frame.create_settings_objs()
    return frame.run_experiments()

if __name__ == "__main__":
    main('parameters.json')
