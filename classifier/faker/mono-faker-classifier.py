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

#     ********************************************************************************************************************************** 
#     *********************** NOTE THAT THIS USES FAKE DATA FOR TESTING PURPOSES ONLY - DO NOT USE ON OUR ACTUAL DATA ******************
#                             This runs on data generated by mono-faker.py as a comparison to our monotherapy simulated
#                             data running through mono-tumor-classifier.py. The tests, as of 6/15/17, came back successfully.
#                             This is simply meant to check our algorithm to ensure that it is running correctly (i.e. learning).
#                             Also note that the dimensions of the arrays are different b/c the fake data is (slightly) different.
#     ********************************************************************************************************************************** 

# How to run
# start tensorflow virtualenv - in OS X:
#   $source [/path/to/tensorflow/bin/activate]
# run mono-faker.py supplying the appropriate args at command line (see mono-faker.py)
# run faker-classifier.py
#   python3 faker-classifier.py

# in mnist.py, mnist.train.images is a tensor (an n-dimensional array)
# with a shape of [55000, 784], where the first dimension is an index
# into the list of images and the second dimension is the index for
# each pixel in each image. in our model, the analogous variable to
# mnist.train.images is a tensor with shape [n, 6] where n is the
# number of data lines in mrhlog.tsv. Therefore, the first dimension
# is an index to the response of a tumor type to a given drug (i.e.
# delta k). The second dimension represents a 11 dimensional array of
# tumor type, which is represented in genetics.tsv by a binary number
# that is then doubled.

# Additionally, in mnist.py, mnist.train.labels is a [55000, 10] array
# of floats, where the second dimension represents a 10 dimensional
# array of labels with one value set to 1 to represent the label for
# the image indexed by the value of the first dimension. In our model,
# the analogous tensor is shape [n, 10] where the second dimension is
# an 11 dimensional array of delta ks (normalized from [-20, 20] to
# [-1, 1])

from os import sys
import tensorflow as tf
import random, time, datetime

global log
path_to_data_file = "inputs/mono-faker-corr-0-05" # SET THIS TO WHATEVER YOU NAMED YOUR FAKE DATA FILE 

#### PARAMETERS ####

class settings:

    def param(self, key):
        return self.params[key]

    params = {} 

    alphabet = {'a0': 0, 'b1': 1, 'c2': 2, 'd3': 3, 'e4': 4, 'f5': 5, 
                 'g6': 6, 'h7': 7, 'i8': 8, 'j9': 9}

    #note these params are for testing fake data only, will give weird results for actual data
    param_specs = {"run_label": ['testing fake data with corr=0.05'],
                   "training_mode": ['normal', 'random'], # normal, stretch, random
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
                   "test_fraction" : [0.1], # DEFAULT : 20 (was 10 in MNIST)
                   "delta_k_training_inclusion_threshold" : [-1], # DEFAULT : -20
                   "delta_k_testing_inclusion_threshold" : [-1], # DEFAULT : 15
                   "normal_training_val": [1]
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
    alpha = split_line[0]
    label = int(split_line[1])
    is_correct_label = bool(int(split_line[2]))

    #input is the alphabet character...
    x = [0 for i in range(10)]
    x[settings.alphabet[alpha]] = 1

    if (settings.param("background_mode") == 'random'):
        y_ = [((random.random()*settings.param("random_background_multiplier"))-settings.param("random_background_shift")) for i in range(10)]
    elif (settings.param("background_mode") == 'constant'):
        y_ = [settings.param("background_constant") for i in range(10)]
    else:
        raise Exception("Bad background_mode")
    #label is label (just a number)
    if (settings.param('training_mode') == 'constant'):
        y_[label] = settings.param('training_constant')
    elif (settings.param('training_mode') == 'random'):
        y_[label] = (random.random()*settings.param('random_training_multiplier'))-settings.param('random_training_shift')
    elif (settings.param('training_mode') == 'normal'):
        y_[label] = settings.param('normal_training_val')
    elif (settings.param('training_mode') == 'stretch'):
        y_[label] = (settings.param('normal_training_val')+settings.param('stretch_training_shift'))/settings.param('stretch_training_divisor')
    else:
        raise Exception("Bad training_mode")

    return x, y_, is_correct_label

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
        x, y_, is_correct_label = parse_line(line, settings)
        training_x.append(x)
        training_y.append(y_)
    print('len training x: ' + str(len(training_x)))
    for line in testing_lines:
        x, y_, is_correct_label = parse_line(line, settings) 
        if is_correct_label:
            testing_x.append(x)
            testing_y.append(y_)
    print('len testing x: ' + str(len(testing_x)))
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
    log = open("runlog-faker.dbl", "a")
    settings.params['ts'] = str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
    print(str(settings.params))
    meta = [] # holds the accuracy from each run
    training_x, training_y, testing_x, testing_y = parse_log(path_to_data_file, settings)

    # Initializing the weights randomly makes no major difference, so long as the stddev is 1.0
    for k in range(settings.param("trial_n")):
        x = tf.placeholder(tf.float32, [None, 10]) #input (10 dim alphabet array)
        W = tf.Variable(tf.zeros([10, 10])) #weights (10 dim alphabet array to 10 dim label array)
        b = tf.Variable(tf.zeros([10])) #biases (10 dim label array)
        
        ###### MAIN DIFFERENCE B/W mnist-adapted.py and mnist-adapted-corr.py ########
        y_ = tf.placeholder(tf.float32, [None, 10]) #"correct" answers (treatments)
        y = tf.matmul(x, W) + b #output (10 wide alphabet array)
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
            batch_xs = [training_x[j] for j in sampled]
            batch_ys = [training_y[j] for j in sampled]
            sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})
        
        # Test 
        correct_prediction = tf.equal(tf.argmax(y, 1), tf.argmax(y_, 1))
        accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
        result = sess.run(accuracy, feed_dict={x: testing_x, y_: testing_y})

        # Record
        meta.append(result)
        settings.params['results'] = meta
        print(str(settings.params))

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
