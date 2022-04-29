# ***** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION *****
# (See ***** below)

#     ************************************************************************************************************************
#     ******************************** REMEMBER TO CHANGE THE EXPERIMENT_LABEL         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************

#     ************************************************************************************************************************
#     ******************************** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************#     ************************************************************************************************************************
#     ******************************** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     ************************************************************************************************************************#     ************************************************************************************************************************
#     ******************************** NOTE THE INITIALIZATION LINES NEEDS TO BE CHANGED FOR PYTHON VERSION    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
import random, time, datetime

global log
path_to_data_file = "inputs/3706274111-mrhlog.xls"

#### PARAMETERS ####

class settings:

    def param(self, key):
        return self.params[key]

    params = {} 

    drug_dict = {'BEVACIZUMAB': 0, 'TEMOZOLOMIDE': 1, 'CABAZITAXEL': 2, 'TCAR': 3, 'DISATINIB': 4, 'NIVOLUMAB': 5, 
                 'DOXORUBICIN': 6, 'DURVALUMAB': 7, 'PEMBROLIZUMAB': 8, 'VARLILUMAB': 9, 'SORAFENIB': 10}

    param_specs = {"run_label": ['testing random weights stddev = 1'],
                   "training_mode": ['normal'], # normal, stretch, random
                   "normal_training_divisor": [20], # 20
                   "stretch_training_shift" : [20.0],
                   "stretch_training_divisor" : [40.0],
                   "training_constant" : [1.0], # training_mode : 'constant' 
                   "random_training_multiplier" : [1.0], # training_mode : 'random'
                   "random_training_shift" : [-2.0], # training_mode : 'random'
                   "background_mode" : ['constant'], # constant v. random
                   "background_constant" : [0.0], # Usually 0.0
                   "random_background_multiplier" : [0.25],
                   "random_background_shift" : [0.125],
                   "GD_step_size" : [0.5], # DEFAULT : 0.5
                   "trial_n" : [5],
                   "test_fraction" : [.20], # DEFAULT : 20 (was 10 in MNIST)
                   "delta_k_training_inclusion_threshold" : [-20], # DEFAULT : -20
                   "delta_k_testing_inclusion_threshold" : [4] # DEFAULT : 15
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

def parse_line(line, settings):
    split_line = line.split('\t')
    tumor_id = int(split_line[1])
    #turn tumor_id into binary representation
    tumor_rep = list("{0:06b}".format(tumor_id//2))
    #tumor_rep is list of strings, convert to int
    for i in range(len(tumor_rep)):
        tumor_rep[i] = int(tumor_rep[i])

    if (settings.param("background_mode") == 'random'):
        y_ = [((random.random()*settings.param("random_background_multiplier"))-settings.param("random_background_shift")) for i in range(11)]
    elif (settings.param("background_mode") == 'constant'):
        y_ = [settings.param("background_constant") for i in range(11)]
    else:
        raise Exception("Bad background_mode")

    drug_name = split_line[3]
    delta_k = int(split_line[-1])

    if (settings.param("training_mode") == 'constant'):
        y_[settings.drug_dict[drug_name]] = settings.param("training_constant")
    elif (settings.param("training_mode") == 'random'):
        y_[settings.drug_dict[drug_name]] = ((random.random()*settings.param("random_training_multiplier"))-settings.param("random_training_shift"))
    elif (settings.param("training_mode") == 'normal'):
        # normalized delta k -- Putatively correct version with 
        y_[settings.drug_dict[drug_name]] = delta_k/settings.param("normal_training_divisor")
    elif (settings.param("training_mode") == 'stretch'):
        y_[settings.drug_dict[drug_name]] = (delta_k+settings.param("stretch_training_shift"))/settings.param("stretch_training_divisor")
    else:
        raise Exception("Bad training_mode")

    return tumor_rep, y_, delta_k

def parse_log(filename, settings):
    training_x = []
    testing_x = []
    training_y = []
    testing_y = []
    dfile = open(filename, 'r')
    lines = dfile.readlines()
    training_lines = lines[int(len(lines)*settings.param("test_fraction")):]
    testing_lines = lines[:int(len(lines)*settings.param("test_fraction"))]
    for line in training_lines:
        x, y_, delta_k = parse_line(line, settings)
        if (delta_k > settings.param("delta_k_training_inclusion_threshold")):
            training_x.append(x)
            training_y.append(y_)
    for line in testing_lines:
        x, y_, delta_k = parse_line(line, settings) 
        if (delta_k > settings.param("delta_k_testing_inclusion_threshold")):
            testing_x.append(x)
            testing_y.append(y_)
    if (len(testing_x) == 0):
        raise Exception("Too few testing samples to run!")

    return training_x, training_y, testing_x, testing_y

#### MAIN FNS ####

#From mnist.py -- adapt this to work with our data from mrhlog.tsv
#
# for _ in range(1000):
#   batch_xs, batch_ys = mnist.train.next_batch(100)
#   sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})
#
#Note that mnist.py uses mnist.train.next_batch to pull up new values for x and
#y_ (correct values for y) - so, need to write small script to parse mrhlog.tsv
#into a format that can be broken down and shoved into the feed_dict from third
#line of above block. Should be fairly straightforward.

def run(settings):
    global log
    log = open("runlog.dbl", "a")
    settings.params['ts'] = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
    print(str(settings.params))
    meta = [] # holds the accuracy from each run
    training_x, training_y, testing_x, testing_y = parse_log(path_to_data_file, settings)

    # Initializing the weights randomly makes no major difference, so long as the stddev is 1.0
    for k in range(settings.param("trial_n")):
        x = tf.placeholder(tf.float32, [None, 6]) #input (tumors)
        W = tf.Variable(tf.truncated_normal([6, 11], mean=0.0, stddev=1.0, dtype=tf.float32, seed=None, name=None)) #weights (6 wide tumor array to 11 wide drug array)
        #W = tf.Variable(tf.zeros([6, 11])) #weights (6 wide tumor array to 11 wide drug array)
        b = tf.Variable(tf.zeros([11])) #biases (11 wide drug array)
        
        ###### MAIN DIFFERENCE B/W mnist-adapted.py and mnist-adapted-corr.py ########
        y_ = tf.placeholder(tf.float32, [None, 11]) #"correct" answers (treatments)
        y = tf.matmul(x, W) + b #output (11 wide drug array)
        cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y))
        ##############################################################################
        
        train_step = tf.train.GradientDescentOptimizer(settings.param("GD_step_size")).minimize(cross_entropy)

        sess = tf.InteractiveSession()
        
        # ***** This is the only thing that doesn't work in both pythons *****
        # There may be some way to conditionalize on this, but for now you just have to edit it.
        tf.global_variables_initializer().run() # Python3
        #tf.initialize_all_variables().run()      # Python 2(.7)

        # Train
        for i in range(1000):
            sampled = random.sample(range(1, len(training_x)), 100)
            batch_xs = [training_x[i] for i in sampled]
            batch_ys = [training_y[i] for i in sampled]
            sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})
        
        # Test 
        correct_prediction = tf.equal(tf.argmax(y, 1), tf.argmax(y_, 1))
        accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
        result = sess.run(accuracy, feed_dict={x: testing_x, y_: testing_y})

        # Record
        meta.append(result)
        settings.params['results'] = meta
        print(str(settings.params))
        log.write(str(settings.params))
        log.write("\n")

    # Final recording
    settings.params['mean_results'] = sum(settings.params['results'])/len(settings.params['results'])
    log.write(str(settings.params))
    log.write("\n")
    log.close()

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

config_and_test()
