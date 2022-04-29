# Implementation of a simple MLP network with one hidden layer. Tested on the iris data set.
# Requires: numpy, sklearn>=0.18.1, tensorflow>=1.0

# NOTE: In order to make the code simple, we rewrite x * W_1 + b_1 = x' * W_1'
# where x' = [x | 1] and W_1' is the matrix W_1 appended with a new row with elements b_1's.
# Similarly, for h * W_2 + b_2
import tensorflow as tf
import numpy as np
from sklearn.model_selection import train_test_split
import random
import json

RANDOM_SEED = random.randint(0, 1000000)
tf.set_random_seed(RANDOM_SEED)

def init_weights(shape):
    """ Weight initialization """
    weights = tf.random_normal(shape, stddev=0.1)
    return tf.Variable(weights)

def forwardprop(X, w_1, w_2):
    """
    Forward-propagation.
    IMPORTANT: yhat is not softmax since TensorFlow's softmax_cross_entropy_with_logits() does that internally.
    """
    h    = tf.nn.sigmoid(tf.matmul(X, w_1))  # The \sigma function
    yhat = tf.matmul(h, w_2)  # The \varphi function
    return yhat

def parse_chunk(chunk):
    minisim_info = json.load(open('max-conv-minisim-information.json', 'r'))
    # minisim_info = json.load(open('minisim-information.json', 'r'))

    # these two values should be the same in all given lines...
    tumor_id = int(chunk[0][2])
    delta_k = float(chunk[0][-1])

    drug_names = [line[4] for line in chunk]

    markers = minisim_info['markers']

    #format the tumor type as a binary number/array for input
    format_string = "{0:0" + str(len(markers)) + "b}"
    x = list(format_string.format(tumor_id))
    #convert list of chars to list of ints
    for i in range(len(x)):
        x[i] = float(x[i])
    x = [1] + x #bias val...


    drug_list = minisim_info['drugs']
    drug_dict = {drug_list[i]: i for i in range(len(drug_list))}


    # for drug_name in drug_names:
    #     y_[drug_dict[drug_name]] = delta_k

    # if delta_k > 0:
    #     for drug_name in drug_names:
    #         y_[drug_dict[drug_name]] = 1
    # elif delta_k < 0:
    #     for drug_name in drug_names:
    #         y_[drug_dict[drug_name]] = -1

    y_ = [0 for i in range(len(drug_dict))]

    for drug in drug_names:
        y_[drug_dict[drug]] = delta_k

    # if drug_names[0] == 'BEVACIZUMAB':
    #     y_ = [delta_k, 0]
    # elif drug_names[0] == 'TEMOZOLOMIDE':
    #     y_ = [0, delta_k]

    return x, y_

def get_mrhlog_data():
    data_file = open('max-conv-minisim-mrhlog.xls', 'r')
    # data_file = open('minisim-mrhlog.xls', 'r')
    lines = data_file.readlines()
    split_lines = [line.split() for line in lines]

    # data_chunks = [[line] for line in split_lines]

    # index 3 is the combo number...
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
    adjusted_std_delta_k = std_delta_k * 1
    print('adjusted-std-delta-k = ',str(adjusted_std_delta_k))
    max_delta_k = max(all_dk)
    min_delta_k = min(all_dk)
    # !!!!!!!!!!! Trying an experiment forcing this to zero mean:
    # mean_delta_k = 0
    dk_positive_thresh = mean_delta_k + adjusted_std_delta_k
    dk_negative_thresh = mean_delta_k - adjusted_std_delta_k
    print(dk_negative_thresh)
    print(dk_positive_thresh)

    # good_chunks = []
    # for chunk in data_chunks:
    #     if int(chunk[0][-1]) >= dk_positive_thresh or int(chunk[0][-1]) <= dk_negative_thresh:
    #         good_chunks.append(chunk)

    good_chunks = data_chunks

    xs = []
    ys = []

    for chunk in good_chunks:
        x, y_ = parse_chunk(chunk)
        xs.append(x)
        ys.append(y_)

    print('len(xs)', len(xs))

    return train_test_split(xs, ys, test_size = 0.2, random_state=RANDOM_SEED)

def main():
    train_X, test_X, train_y, test_y = get_mrhlog_data()

    # Layer's sizes
    x_size = len(train_X[0])   # Number of input nodes: 1 bias, 6 mutation features
    h_size = 50             # Number of hidden nodes
    y_size = len(train_y[0])   # Number of outcomes (11 drugs)

    new_order_training = random.sample(range(len(train_X)), len(train_X))
    new_order_testing = random.sample(range(len(test_X)), len(test_X))

    train_X = [train_X[i] for i in new_order_training]
    train_y = [train_y[i] for i in new_order_training]
    test_X = [test_X[i] for i in new_order_testing]
    test_y = [test_y[i] for i in new_order_testing]

    print('x_size', x_size)
    print('y_size', y_size)

    # Symbols
    X = tf.placeholder("float", shape=[None, x_size])
    y = tf.placeholder("float", shape=[None, y_size])

    # Weight initializations
    w_1 = init_weights((x_size, h_size))
    w_2 = init_weights((h_size, y_size))

    # Forward propagation
    yhat    = forwardprop(X, w_1, w_2)
    predict = tf.argmax(yhat, axis=1)

    #trying to do double drug accuracy...
    # predict = tf.nn.top_k(yhat, 2)

    # max_yhat = tf.nn.top_k(yhat, 2).indices
    # max_max_yhat = tf.nn.top_k(max_yhat, 2).values

    # max_y_train = tf.nn.top_k(train_y, 2).indices
    # max_max_y_train = tf.nn.top_k(max_y_train, 2).values

    # max_y_test = tf.nn.top_k(test_y, 2).indices
    # max_max_y_test = tf.nn.top_k(max_y_test, 2).values

    # Backward propagation
    cost    = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y, logits=yhat))
    updates = tf.train.GradientDescentOptimizer(0.01).minimize(cost)

    # Run SGD
    sess = tf.Session()
    init = tf.global_variables_initializer()
    sess.run(init)

    for epoch in range(30):
        # Train with each example

        for i in range(len(train_X)):
            sess.run(updates, feed_dict={X: train_X[i: i + 1], y: train_y[i: i + 1]})

        # train_accuracy = np.mean(max_max_y_train == sess.run(max_max_yhat, feed_dict={X: train_X, y: train_y}))
        # test_accuracy = np.mean(max_max_y_test == sess.run(max_max_yhat, feed_dict={X: train_X, y: train_y}))

        # train_accuracy = np.mean( == sess.run(predict, feed_dict={X: train_X, y: train_y}))
        # test_accuracy = np.mean(max_y_test == sess.run(predict, feed_dict={X: test_X, y: test_y}))

        train_accuracy = np.mean(np.argmax(train_y, axis=1) ==
                                 sess.run(predict, feed_dict={X: train_X, y: train_y}))
        test_accuracy  = np.mean(np.argmax(test_y, axis=1) ==
                                 sess.run(predict, feed_dict={X: test_X, y: test_y}))

        print("Epoch = %d, train accuracy = %.2f%%, test accuracy = %.2f%%"
              % (epoch + 1, 100. * train_accuracy, 100. * test_accuracy))

    sess.close()

if __name__ == '__main__':
    main()