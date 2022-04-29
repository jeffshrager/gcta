# ***** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION *****
# (See ***** below)

#     ************************************************************************************************************************
#     ******************************** REMEMBER TO CHANGE THE EXPERIMENT_LABEL         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************

#     ************************************************************************************************************************
#     ******************************** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************#    
#     ******************************** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************#   
#     ******************************** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************


#     ************************************************************************************************************************
#     ******************************** CHANGE THE DIVIDE BY 2 ON THE TUMOR CODE ONCE SIM IS UPDATED     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************# 
#     ******************************** CHANGE THE DIVIDE BY 2 ON THE TUMOR CODE ONCE SIM IS UPDATED   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************#    
#     ******************************** CHANGE THE DIVIDE BY 2 ON THE TUMOR CODE ONCE SIM IS UPDATED    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

# in mnist.py, mnist.train.images is a tensor (an n-dimensional array)
# with a shape of [55000, 784], where the first dimension is an index
# into the list of images and the second dimension is the index for
# each pixel in each image. in our model, the analogous variable to
# mnist.train.images is a tensor with shape [n, 6] where n is the
# number of data lines in mrhlog.tsv. Therefore, the first dimension
# is an index to the response of a tumor type to a given drug (i.e.
# delta k). The second dimension represents a 6 dimensional array of
# tumor type, which is represented in genetics.tsv by a binary number
# that is then doubled.

# Additionally, in mnist.py, mnist.train.labels is a [55000, 10] array
# of floats, where the second dimension represents a 10 dimensional
# array of labels with one value set to 1 to represent the label for
# the image indexed by the value of the first dimension. In our model,
# the analogous tensor is shape [n, 11] where the second dimension is
# an 11 dimensional array of delta ks (normalized from [-20, 20] to
# [-1, 1])

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

    #bevacizumab 96 (48 = 110000) EGFR IDH1

    drug_dict = {'BEVACIZUMAB': 0, 'TEMOZOLOMIDE': 1, 'CABAZITAXEL': 2, 'TCAR': 3, 'DISATINIB': 4, 'NIVOLUMAB': 5, 
                 'DOXORUBICIN': 6, 'DURVALUMAB': 7, 'PEMBROLIZUMAB': 8, 'VARLILUMAB': 9, 'SORAFENIB': 10}

    drug_combo_list_dict = {}

    drug_combos = [dc for dc in combos(drug_dict.keys(), 2)]
    for i in range(len(drug_combos)):
        drug_combo_list_dict[i] = list(drug_combos[i])

    for i in range(len(drug_dict.keys())):
        drug_combo_list_dict[i+len(drug_combos)] = [list(drug_dict.keys())[i], list(drug_dict.keys())[i]]

    for i in drug_combo_list_dict.keys():
        drug_combo_list_dict[i].sort()

    # print(drug_combo_list_dict)

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
                   "trial_n" : [3], # DEFAULT : 1
                   "test_fraction" : [0.1], # DEFAULT : 0.2
                   "delta_k_training_inclusion_threshold" : [10], # DEFAULT : -20
                   "delta_k_testing_inclusion_threshold" : [10], # DEFAULT : -20
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

    tumor_id = int(split_line_one[2])
    #turn tumor_id into binary representation
    x = list("{0:06b}".format(tumor_id//2))
    #x is list of strings, convert to int
    for i in range(len(x)):
        x[i] = int(x[i])

    if (settings.param("background_mode") == 'random'):
        y_ = [((random.random()*settings.param("random_background_multiplier"))-settings.param("random_background_shift")) for i in range(66)]
    elif (settings.param("background_mode") == 'constant'):
        y_ = [settings.param("background_constant") for i in range(66)]
    else:
        raise Exception("Bad background_mode")

    drug_name_one = split_line_one[4]
    drug_name_two = split_line_two[4]
    drug_name_list = [drug_name_one, drug_name_two]
    drug_name_list.sort()

    index_of_drug_cocktail = None
    for i in range(len(settings.drug_combo_list_dict)):
        if settings.drug_combo_list_dict[i] == drug_name_list:
            index_of_drug_cocktail = i
            break

    #note delta_k should be same for line one and line two (that's the point)
    delta_k = int(split_line_one[-1])

    if (settings.param("training_mode") == 'constant'):
        if delta_k > 0:
            y_[index_of_drug_cocktail] = settings.param("training_constant_positive")
        elif delta_k < 0: 
            y_[index_of_drug_cocktail] = settings.param("training_constant_negative")
        else:
            y_[index_of_drug_cocktail] = 0
    elif (settings.param("training_mode") == 'random'):
        y_[index_of_drug_cocktail] = ((random.random()*settings.param("random_training_multiplier"))-settings.param("random_training_shift"))    
    elif (settings.param("training_mode") == 'true_random'):
        y_ = [(random.random()*settings.param('random_training_multiplier'))-settings.param('random_training_shift') for i in range(66)]
    elif (settings.param("training_mode") == 'normal'):
        # normalized delta k -- Putatively correct version with 
        y_[index_of_drug_cocktail] = delta_k
    elif (settings.param("training_mode") == 'stretch'):
        y_[index_of_drug_cocktail] = (delta_k+settings.param("stretch_training_shift"))/settings.param("stretch_training_divisor")
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

def run(settings):
    global log
    settings.params['ts'] = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
    print(str(settings.params))
    meta = [] # holds the accuracy from each run
    training_x, training_y, testing_x, testing_y = parse_log(path_to_data_file, settings)

    for k in range(settings.param("trial_n")):
        x = tf.placeholder(tf.float32, [None, 6]) #input (tumors)
        W = tf.Variable(tf.zeros([6, 66]))
        # W = generate_weight_matrix(settings)
        b = tf.Variable(tf.zeros([66]))
    
        ###### MAIN DIFFERENCE B/W mnist-adapted.py and mnist-adapted-corr.py ########
        y_ = tf.placeholder(tf.float32, [None, 66]) #"correct" answers (treatments)
        y = tf.matmul(x, W) + b #output (11 or 22 wide drug array)
        cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y))
        ##############################################################################
        
        train_step = tf.train.GradientDescentOptimizer(settings.param("GD_step_size")).minimize(cross_entropy)

        sess = tf.InteractiveSession()
        
        # ***** This is the only thing that doesn't work in both pythons *****
        # There may be some way to conditionalize on this, but for now you just have to edit it.
        tf.global_variables_initializer().run() # Python3
        #tf.initialize_all_variables().run()      # Python 2(.7)

        # Train
        for i in range(settings.param("num_training_cycles")):
            sampled = random.sample(range(1, len(training_x)), settings.param("num_training_samples"))
            batch_xs = [training_x[i] for i in sampled]
            batch_ys = [training_y[i] for i in sampled]
            sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})

        # Test 
        max_y = tf.nn.top_k(y, 2).indices #what the network thinks is going on
        max_y_ = tf.nn.top_k(y_, 2).indices #what's really going on -- target
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
        run(settings())
        # try:
        #     run(settings())
        # except Exception as exc:
        #     pnl(str(exc))
        #     pnl("Failed! Moving on...")

#### EXECUTION ####

log = open("runlog.dbl", "a")

config_and_test()

log.close()
