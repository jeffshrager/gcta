from os import sys
import tensorflow as tf
import numpy as np
import random
import datetime
import time
import json
from itertools import combinations as combos, product as iterproduct

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
        mut_dict = dict((v,k) for k,v in self.settings.markers_dict.items())
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
        self.settings.params['testing-length'] = len(testing_x)
        self.settings.params['training-length'] = len(training_x)

        meta_results = []
        for k in range(self.settings.param("trial-n")):
            g = tf.Graph()

            with g.as_default():
                tf.set_random_seed(random.randint(0, 10000000))

                sess = tf.InteractiveSession()    

                x = tf.placeholder(tf.float32, [None, self.num_markers]) #input (tumors)
                W = self.generate_tf_variable([self.num_markers, self.num_drugs], 'noisy')
                      
                b = self.generate_tf_variable([self.num_drugs], 'zeros')

                y = tf.matmul(x, W) + b

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
                    sampled = random.sample(range(1, len(training_x)), self.settings.param("num-training-samples"))
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
                print('On trial', str(k+1) + '/' + str(self.settings.param('trial-n')), 'of experiment configuration', str(settings_obj_num) + '.')

        self.settings.params['results'] = meta_results
        self.settings.params['mean-results'] = sum(meta_results)/len(meta_results)
        self.settings.hypotheses = self.generate_hypotheses()

        return self.settings.params, self.settings.hypotheses

    def parse_log(self):
        training_x = []
        testing_x = []
        training_y = []
        testing_y = []

        self.data_file.seek(0)
        lines = self.data_file.readlines()

        training_lines = lines[:int(len(lines)*(1-self.settings.param("test-fraction")))]
        testing_lines = lines[int(len(lines)*(1-self.settings.param("test-fraction"))):]

        #loop through training lines to produce training_x, training_y
        i = 0
        while i < len(training_lines)-1:
            x, y_, delta_k = self.parse_two_lines(training_lines[i], training_lines[i+1])
            if abs(delta_k) >= self.settings.param("delta-k-training-inclusion-abs-threshold"):
                training_x.append(x)
                training_y.append(y_)
            i += 2

        #loop through testing lines to produce testing_x, testing_y
        i = 0
        while i < len(testing_lines)-1:
            x, y_, delta_k = self.parse_two_lines(testing_lines[i], testing_lines[i+1])
            if delta_k >= self.settings.param("delta-k-testing-inclusion-abs-threshold"):
                testing_x.append(x)
                testing_y.append(y_)
            i += 2

        if (len(testing_x) == 0):
            raise Exception("Too few testing samples to run!")

        return training_x, training_y, testing_x, testing_y

    def parse_two_lines(self, line_one, line_two):
        split_line_one = line_one.split()
        split_line_two = line_two.split()

        #these two values should be the same for both lines
        tumor_id = int(split_line_one[2])
        delta_k = int(split_line_one[-1])

        drug_name_one = split_line_one[4]
        drug_name_two = split_line_two[4]

        #turn tumor_id into binary representation
        x = list("{0:06b}".format(tumor_id))
        #x is list of strings, convert to int
        for i in range(len(x)):
            x[i] = int(x[i])

        if self.settings.param("background-mode") == 'random':
            y_ = [((random.random()*self.settings.param("random-background-multiplier"))-self.settings.param("random-background-shift")) for i in range(self.num_drugs)]
        elif self.settings.param("background-mode") == 'constant':
            y_ = [self.settings.param("background-constant") for i in range(self.num_drugs)]
        else:
            raise Exception("Bad background-mode :(")

        if (self.settings.param("training-mode") == 'constant'):
            if delta_k > 0:
                y_[self.settings.drug_dict[drug_name_one]] = self.settings.param("training-constant-positive")
                y_[self.settings.drug_dict[drug_name_two]] = self.settings.param("training-constant-positive")
            elif delta_k < 0: 
                y_[self.settings.drug_dict[drug_name_one]] = self.settings.param("training-constant-negative")
                y_[self.settings.drug_dict[drug_name_two]] = self.settings.param("training-constant-negative")
            else:
                y_[self.settings.drug_dict[drug_name_one]] = 0
                y_[self.settings.drug_dict[drug_name_two]] = 0
        elif (self.settings.param("training-mode") == 'random'):
            y_[self.settings.drug_dict[drug_name_one]] = ((random.random()*self.settings.param("random-training-multiplier"))-self.settings.param("random-training-shift"))
            y_[self.settings.drug_dict[drug_name_two]] = ((random.random()*self.settings.param("random-training-multiplier"))-self.settings.param("random-training-shift"))    
        elif (self.settings.param("training-mode") == 'true-random'):
            #note that this erases whatever default values are set with self.settings.param("background_mode")... although it is truly random noise for all values
            y_ = [(random.random()*self.settings.param('random-training-multiplier'))-self.settings.param('random-training-shift') for i in range(self.num_drugs)]
        elif (self.settings.param("training-mode") == 'normal'):
            # "normalized" delta k (i.e. divide by delta_k scale size)
            y_[self.settings.drug_dict[drug_name_one]] = delta_k/self.settings.param("normal-training-divisor")
            y_[self.settings.drug_dict[drug_name_two]] = delta_k/self.settings.param("normal-training-divisor")
        elif (self.settings.param("training-mode") == 'stretch'):
            y_[self.settings.drug_dict[drug_name_one]] = (delta_k+self.settings.param("stretch-training-shift"))/self.settings.param("stretch-training-divisor")     
            y_[self.settings.drug_dict[drug_name_two]] = (delta_k+self.settings.param("stretch-training-shift"))/self.settings.param("stretch-training-divisor")
        else:
            raise Exception("Bad training-mode")

        return x, y_, delta_k

class ExpSettings:
    def __init__(self, keys, values, drug_dict, markers_dict):
        self.params = dict(zip(keys, values))
        self.drug_dict = drug_dict
        self.markers_dict = markers_dict
        self.params['best-accuracy'] = 0
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
        self.run_log = open('experiments/'+self.run_label+'/'+self.log_name, 'a')
        self.path_to_data_file = 'experiments/'+self.run_label+'/'+self.data_name
        print('data file: ' + self.path_to_data_file)
        self.data_file = open(self.path_to_data_file)
        self.drugs = list(self.params['simulator']['drugs'].keys())
        self.drugs = [drug.upper() for drug in self.drugs]
        self.drug_dict = dict(zip(self.drugs, range(len(self.drugs))))
        self.markers = self.params['simulator']['markers']

    def pnl(self, to_pnl):
        self.run_log.write(to_pnl+'\n')
        print(to_pnl)

    def create_settings_objs(self):
        keys = list(self.learner_params.keys())
        raw_vals = list(self.learner_params.values())
        params_products = iterproduct(*raw_vals)

        for prod in params_products:
            vals = [v for v in prod]
            self.settings_objs.append(ExpSettings(keys, vals, self.drug_dict, self.markers))

    def run_experiments(self):
        for i in range(len(self.settings_objs)):
            settings_obj = self.settings_objs[i]
            l = Learner(settings_obj, self.data_file)
            to_log, hypotheses = l.run(i+1)
            to_log['hypotheses'] = str(hypotheses).replace('\'', '"').upper()
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