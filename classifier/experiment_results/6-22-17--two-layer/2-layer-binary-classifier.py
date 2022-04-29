# ***** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION *****
# (See ***** below)

#     ************************************************************************************************************************
#     ******************************** REMEMBER TO CHANGE THE EXPERIMENT_LABEL         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************


# How to run
# start tensorflow virtualenv - in OS X:
#   $source [/path/to/tensorflow/bin/activate]
# run mnist-adapted-corr.py, provide path to mrhlog.tsv file with input data
#   python3 mnist-adapted.py [/path/to/mrhlog.tsv]

# Alt (Jeff's preferred path): 
#   $source ~/tensorflow/bin/activate
#   $python (or python3, note having to change for py version)
#   >>> execfile("mnist-adapted-corr.py")

from os import sys
import tensorflow as tf
import numpy as np
import random, time, datetime
from itertools import combinations as combos

global log
path_to_data_file = "inputs/3706550898-mrhlog.xls"

#### PARAMETERS ####

class settings:

    def param(self, key):
        return self.params[key]

    params = {} 

    #drug codes like so : bevacizumab 96 (48 = 110000) EGFR IDH1

    drug_dict = {'BEVACIZUMAB': 0, 'TEMOZOLOMIDE': 1, 'CABAZITAXEL': 2, 'TCAR': 3, 'DISATINIB': 4, 'NIVOLUMAB': 5, 
                 'DOXORUBICIN': 6, 'DURVALUMAB': 7, 'PEMBROLIZUMAB': 8, 'VARLILUMAB': 9, 'SORAFENIB': 10}

    param_specs = {"run_label": ['test'],
                   "training_mode": ['constant', 'true_random'], # constant, normal, stretch, random, true_random
                   "normal_training_divisor": [20.0], # 20
                   "stretch_training_shift" : [20.0],
                   'bias_val' : [0],
                   "stretch_training_divisor"  : [40.0],
                   "training_constant_positive" : [1.0], # training_mode : 'constant'
                   "training_constant_negative" : [-1.0], # training_mode : 'constant' 
                   "random_training_multiplier" : [2.0], # training_mode : 'random'
                   "random_training_shift" : [1.0], # training_mode : 'random'
                   "background_mode" : ['constant'], # constant v. random
                   "background_constant" : [0.0], # Usually 0.0
                   "random_background_multiplier" : [0.5], #.25 before
                   "random_background_shift" : [0.25], #.125 before
                   "GD_step_size" : [0.3], # DEFAULT : 0.5
                   "trial_n" : [5], # DEFAULT : 1
                   "test_fraction" : [0.1], # DEFAULT : 0.2
                   "delta_k_training_inclusion_abs_threshold" : [15], # DEFAULT : -20
                   "delta_k_testing_inclusion_abs_threshold" : [15], # DEFAULT : -20
                   "hypothesis" : ['all_zeros'], #'all_zeros', 'random_hypothesis', 'bias_hypothesis'
                   "num_training_cycles": [1000],
                   "num_training_samples": [100]
                   }

global param_specs_keys
param_specs_keys = list(settings.param_specs.keys())

#### UTILS ####

def pnl(str): # print and log
    global log
    print(str)
    log.write(str)
    log.write('\n')

#### DATA IMPORT AND INPUTS SETUP ####

def parse_two_lines(line_one, line_two, settings):
    split_line_one = line_one.split('\t')
    split_line_two = line_two.split('\t')

    #these two values should be the same for both lines
    assert(split_line_one[-1] == split_line_two[-1])
    assert(split_line_one[2] == split_line_two[2])
    tumor_id = int(split_line_one[2])
    delta_k = int(split_line_one[-1])

    drug_name_one = split_line_one[4]
    drug_name_two = split_line_two[4]

    #turn tumor_id into binary representation
    x = list("{0:06b}".format(tumor_id//2))
    #x is list of strings, convert to int
    for i in range(len(x)):
        x[i] = int(x[i])

    if settings.param("background_mode") == 'random':
        y_ = [((random.random()*settings.param("random_background_multiplier"))-settings.param("random_background_shift")) for i in range(11)]
    elif settings.param("background_mode") == 'constant':
        y_ = [settings.param("background_constant") for i in range(11)]
    else:
        raise Exception("Bad background_mode :(")

    if (settings.param("training_mode") == 'constant'):
        if delta_k > 0:
            y_[settings.drug_dict[drug_name_one]] = settings.param("training_constant_positive")
            y_[settings.drug_dict[drug_name_two]] = settings.param("training_constant_positive")
        elif delta_k < 0: 
            y_[settings.drug_dict[drug_name_one]] = settings.param("training_constant_negative")
            y_[settings.drug_dict[drug_name_two]] = settings.param("training_constant_negative")
        else:
            y_[settings.drug_dict[drug_name_one]] = 0
            y_[settings.drug_dict[drug_name_two]] = 0
    elif (settings.param("training_mode") == 'random'):
        y_[settings.drug_dict[drug_name_one]] = ((random.random()*settings.param("random_training_multiplier"))-settings.param("random_training_shift"))
        y_[settings.drug_dict[drug_name_two]] = ((random.random()*settings.param("random_training_multiplier"))-settings.param("random_training_shift"))    
    elif (settings.param("training_mode") == 'true_random'):
        #note that this erases whatever default values are set with settings.param("background_mode")... although it is truly random noise for all values
        y_ = [(random.random()*settings.param('random_training_multiplier'))-settings.param('random_training_shift') for i in range(11)]
    elif (settings.param("training_mode") == 'normal'):
        # "normalized" delta k (i.e. divide by delta_k scale size)
        y_[settings.drug_dict[drug_name_one]] = delta_k/settings.param("normal_training_divisor")
        y_[settings.drug_dict[drug_name_two]] = delta_k/settings.param("normal_training_divisor")
    elif (settings.param("training_mode") == 'stretch'):
        y_[settings.drug_dict[drug_name_one]] = (delta_k+settings.param("stretch_training_shift"))/settings.param("stretch_training_divisor")     
        y_[settings.drug_dict[drug_name_two]] = (delta_k+settings.param("stretch_training_shift"))/settings.param("stretch_training_divisor")
    else:
        raise Exception("Bad training_mode")

    return x, y_, delta_k

def parse_log(filename, settings):
    training_x = []
    testing_x = []
    training_y = []
    testing_y = []

    dfile = open(filename, 'r')
    lines = dfile.readlines()

    training_lines = lines[int(len(lines)*settings.param("test_fraction")):]
    testing_lines = lines[:int(len(lines)*settings.param("test_fraction"))]

    #loop through training lines to produce training_x, training_y
    i = 0
    while i < len(training_lines)-1:
        x, y_, delta_k = parse_two_lines(training_lines[i], training_lines[i+1], settings)
        if abs(delta_k) >= settings.param("delta_k_training_inclusion_abs_threshold"):
            training_x.append(x)
            training_y.append(y_)
        i += 2

    #loop through testing lines to produce testing_x, testing_y
    i = 0
    while i < len(testing_lines)-1:
        x, y_, delta_k = parse_two_lines(testing_lines[i], testing_lines[i+1], settings)
        if delta_k > settings.param("delta_k_testing_inclusion_abs_threshold"):
            testing_x.append(x)
            testing_y.append(y_)
        i += 2

    if (len(testing_x) == 0):
        raise Exception("Too few testing samples to run!")

    return training_x, training_y, testing_x, testing_y

#### MAIN FNS ####

def generate_tf_variable(shape, mode):
    if mode == 'zeros':
        return tf.Variable(tf.zeros(shape))
    elif mode == 'noisy':
        return tf.Variable(tf.truncated_normal(shape, stddev=0.1))
    else:
        raise Exception("Bad mode for tf variable generation")

def run(settings):
    global log
    settings.params['ts'] = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
    print(str(settings.params))
    meta = [] # holds the accuracy from each run
    training_x, training_y, testing_x, testing_y = parse_log(path_to_data_file, settings)

    for k in range(settings.param("trial_n")):

        g = tf.Graph()

        with g.as_default():
            tf.set_random_seed(random.randint(0, 1000000000))
            
            #two layer network - first layer == two units, second layer just adds and softmaxes
            x = tf.placeholder(tf.float32, [None, 6]) #input (tumors)

            #first layer -- two units
            W_1 = generate_tf_variable([6, 11], 'noisy')
            W_2 = generate_tf_variable([6, 11], 'noisy')

            b_1 = generate_tf_variable([11], 'noisy')
            b_2 = generate_tf_variable([11], 'noisy')

            y_ = tf.placeholder(tf.float32, [None, 11]) #Treatments -- "correct answers"

            y_1 = tf.matmul(x, W_1) + b_1
            y_2 = tf.matmul(x, W_2) + b_2

            #second layer -- sum of first layer's two units
            y = y_1 + y_2

            cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y))
            train_step = tf.train.GradientDescentOptimizer(settings.param("GD_step_size")).minimize(cross_entropy)

            sess = tf.InteractiveSession()
            
            # ***** This is the only thing that doesn't work in both pythons *****
            # NEW: check python version and run according initializer for tf
            if sys.version[0] == '2':
                tf.initialize_all_variables().run()
            else:
                tf.global_variables_initializer().run()

            # Train
            for i in range(settings.param("num_training_cycles")):
                sampled = random.sample(range(1, len(training_x)), settings.param("num_training_samples"))
                batch_xs = [training_x[i] for i in sampled]
                batch_ys = [training_y[i] for i in sampled]
                sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})

            # Test 
            max_y = tf.nn.top_k(y, 2).indices
            max_y_ = tf.nn.top_k(y_, 2).indices
            max_max_y = tf.nn.top_k(max_y, 2).values
            max_max_y_ = tf.nn.top_k(max_y_, 2).values
            correct_prediction = tf.equal(max_max_y, max_max_y_)

            # correct_prediction = tf.equal(tf.nn.top_k(y, 2).indices, tf.nn.top_k(y_, 2).indices)
            accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
            result = sess.run(accuracy, feed_dict={x: testing_x, y_: testing_y})

            # Record
            meta.append(result)
            settings.params['results'] = meta
            print(str(settings.params))

    # Final recording
    settings.params['mean_results'] = sum(settings.params['results'])/len(settings.params['results'])
    settings.params['testing_length'] = len(testing_y) #...
    settings.params['training_length'] = len(training_y)
    log.write(str(settings.params))
    log.write("\n")
    # log.close()

def config_and_test(index=0):
    global param_specs_keys
    if index < len(param_specs_keys): 
        for param_value in settings.param_specs[param_specs_keys[index]]:
            settings.params[param_specs_keys[index]] = param_value
            config_and_test(index + 1)
    else:
        try:
            run(settings())
        except Exception as exc:
            pnl(str(exc))
            pnl("Failed! Moving on...")

#### EXECUTION ####

log = open("runlog.dbl", "a")

config_and_test()

log.close()
