##########################################################
## ccAF:  actinn.py                                     ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @Adapted from:  https://github.com/mafeiyang/ACTINN/ ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

##########################################
## Load Python packages for classifiers ##
##########################################

# General
import math
import numpy as np
import pandas as pd

# Classifier specific
import tensorflow as tf
from tensorflow.python.framework import ops

###############
## Functions ##
###############

# Get common genes, normalize  and scale the sets
def scale_single_set(data_set, train = False, logData = False):
    # input -- set to be scaled
    # output -- scaled sets
    genes = data_set.index
    total_set = np.array(data_set, dtype=np.float32)
    total_set = np.divide(total_set, np.sum(total_set, axis=0, keepdims=True)) * 10000
    if logData:
        total_set = np.log2(total_set+1)
    if train==True:
        expr = np.sum(total_set, axis=1)
        total_set = total_set[np.logical_and(expr >= np.percentile(expr, 1), expr <= np.percentile(expr, 99)),]
        genes = genes[np.logical_and(expr >= np.percentile(expr, 1), expr <= np.percentile(expr, 99))]
        cv = np.std(total_set, axis=1) / np.mean(total_set, axis=1)
        total_set = total_set[np.logical_and(cv >= np.percentile(cv, 1), cv <= np.percentile(cv, 99)),]
        genes = genes[np.logical_and(cv >= np.percentile(cv, 1), cv <= np.percentile(cv, 99))]
        return (total_set, genes)
    else:
        return total_set

# Get common genes, normalize  and scale the sets
def scale_sets(sets, logData = False):
    # input -- a list of all the sets to be scaled
    # output -- scaled sets
    common_genes = set(sets[0].index)
    for i in range(1, len(sets)):
        common_genes = set.intersection(set(sets[i].index),common_genes)
    common_genes = sorted(list(common_genes))
    sep_point = [0]
    for i in range(len(sets)):
        sets[i] = sets[i].loc[common_genes,]
        sep_point.append(sets[i].shape[1])
    total_set = np.array(pd.concat(sets, axis=1, sort=False), dtype=np.float32)
    total_set = np.divide(total_set, np.sum(total_set, axis=0, keepdims=True)) * 10000
    if logData:
        total_set = np.log2(total_set+1)
    expr = np.sum(total_set, axis=1)
    total_set = total_set[np.logical_and(expr >= np.percentile(expr, 1), expr <= np.percentile(expr, 99)),]
    cv = np.std(total_set, axis=1) / np.mean(total_set, axis=1)
    total_set = total_set[np.logical_and(cv >= np.percentile(cv, 1), cv <= np.percentile(cv, 99)),]
    for i in range(len(sets)):
        sets[i] = total_set[:, sum(sep_point[:(i+1)]):sum(sep_point[:(i+2)])]
    return sets

# Turn labels into matrix
def one_hot_matrix(labels, C):
    # input -- labels (true labels of the sets), C (# types)
    # output -- one hot matrix with shape (# types, # samples)
    C = tf.constant(C, name = "C")
    one_hot_matrix = tf.one_hot(labels, C, axis = 0)
    sess = tf.Session()
    one_hot = sess.run(one_hot_matrix)
    sess.close()
    return one_hot

# Function to create placeholders
def create_placeholders(n_x, n_y):
    X = tf.placeholder(tf.float32, shape = (n_x, None))
    Y = tf.placeholder(tf.float32, shape = (n_y, None))
    return X, Y

# Initialize parameters
def initialize_parameters(nf, ln1, ln2, ln3, nt):
    # input -- nf (# of features), ln1 (# nodes in layer1), ln2 (# nodes in layer2), nt (# types)
    # output -- a dictionary of tensors containing W1, b1, W2, b2, W3, b3
    tf.set_random_seed(3) # set seed to make the results consistant
    W1 = tf.get_variable("W1", [ln1, nf], initializer = tf.contrib.layers.xavier_initializer(seed = 3))
    b1 = tf.get_variable("b1", [ln1, 1], initializer = tf.zeros_initializer())
    W2 = tf.get_variable("W2", [ln2, ln1], initializer = tf.contrib.layers.xavier_initializer(seed = 3))
    b2 = tf.get_variable("b2", [ln2, 1], initializer = tf.zeros_initializer())
    W3 = tf.get_variable("W3", [ln3, ln2], initializer = tf.contrib.layers.xavier_initializer(seed = 3))
    b3 = tf.get_variable("b3", [ln3, 1], initializer = tf.zeros_initializer())
    W4 = tf.get_variable("W4", [nt, ln3], initializer = tf.contrib.layers.xavier_initializer(seed = 3))
    b4 = tf.get_variable("b4", [nt, 1], initializer = tf.zeros_initializer())
    parameters = {"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}
    return parameters

# Forward propagation function
def forward_propagation(X, parameters):
    # function model: LINEAR -> RELU -> LINEAR -> RELU -> LINEAR -> SOFTMAX
    # input -- dataset with shape (# features, # sample), parameters "W1", "b1", "W2", "b2", "W3", "b3"
    # output -- the output of the last linear unit
    W1 = parameters['W1']
    b1 = parameters['b1']
    W2 = parameters['W2']
    b2 = parameters['b2']
    W3 = parameters['W3']
    b3 = parameters['b3']
    W4 = parameters['W4']
    b4 = parameters['b4']
    # forward calculations
    Z1 = tf.add(tf.matmul(W1, X), b1)
    A1 = tf.nn.relu(Z1)
    Z2 = tf.add(tf.matmul(W2, A1), b2)
    A2 = tf.nn.relu(Z2)
    Z3 = tf.add(tf.matmul(W3, A2), b3)
    A3 = tf.nn.relu(Z3)
    Z4 = tf.add(tf.matmul(W4, A3), b4)
    return Z4

# Compute cost
def compute_cost(Z4, Y, parameters, lambd=0.01):
    # input -- Z3 (output of forward propagation with shape (# types, # samples)), Y (true labels, same shape as Z3)
    # output -- tensor of teh cost function
    logits = tf.transpose(Z4)
    labels = tf.transpose(Y)
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits = logits, labels = labels)) + \
    (tf.nn.l2_loss(parameters["W1"]) + tf.nn.l2_loss(parameters["W2"]) + tf.nn.l2_loss(parameters["W3"]) + tf.nn.l2_loss(parameters["W4"])) * lambd
    return cost

# Get the mini batches
def random_mini_batches(X, Y, mini_batch_size=32, seed=1):
    # input -- X (training set), Y (true labels)
    # output -- mini batches
    ns = X.shape[1]
    mini_batches = []
    np.random.seed(seed)
    # shuffle (X, Y)
    permutation = list(np.random.permutation(ns))
    shuffled_X = X[:, permutation]
    shuffled_Y = Y[:, permutation]
    # partition (shuffled_X, shuffled_Y), minus the end case.
    num_complete_minibatches = int(math.floor(ns/mini_batch_size)) # number of mini batches of size mini_batch_size in your partitionning
    for k in range(0, num_complete_minibatches):
        mini_batch_X = shuffled_X[:, k * mini_batch_size : k * mini_batch_size + mini_batch_size]
        mini_batch_Y = shuffled_Y[:, k * mini_batch_size : k * mini_batch_size + mini_batch_size]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    # handling the end case (last mini-batch < mini_batch_size)
    if ns % mini_batch_size != 0:
        mini_batch_X = shuffled_X[:, num_complete_minibatches * mini_batch_size : ns]
        mini_batch_Y = shuffled_Y[:, num_complete_minibatches * mini_batch_size : ns]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    return mini_batches

# Forward propagation for prediction
def forward_propagation_for_predict(X, parameters):
    # input -- X (dataset used to make prediction), papameters after training
    # output -- the output of the last linear unit
    W1 = parameters['W1']
    b1 = parameters['b1']
    W2 = parameters['W2']
    b2 = parameters['b2']
    W3 = parameters['W3']
    b3 = parameters['b3']
    W4 = parameters['W4']
    b4 = parameters['b4']
    Z1 = tf.add(tf.matmul(W1, X), b1)
    A1 = tf.nn.relu(Z1)
    Z2 = tf.add(tf.matmul(W2, A1), b2)
    A2 = tf.nn.relu(Z2)
    Z3 = tf.add(tf.matmul(W3, A2), b3)
    A3 = tf.nn.relu(Z3)
    Z4 = tf.add(tf.matmul(W4, A3), b4)
    return Z4

# Predict function
def predict(X, parameters):
    # input -- X (dataset used to make prediction), papameters after training
    # output -- prediction
    W1 = tf.convert_to_tensor(parameters["W1"])
    b1 = tf.convert_to_tensor(parameters["b1"])
    W2 = tf.convert_to_tensor(parameters["W2"])
    b2 = tf.convert_to_tensor(parameters["b2"])
    W3 = tf.convert_to_tensor(parameters["W3"])
    b3 = tf.convert_to_tensor(parameters["b3"])
    W4 = tf.convert_to_tensor(parameters["W4"])
    b4 = tf.convert_to_tensor(parameters["b4"])
    params = {"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}
    x = tf.placeholder("float")
    z4 = forward_propagation_for_predict(x, params)
    p = tf.argmax(z4)
    sess = tf.Session()
    prediction = sess.run(p, feed_dict = {x: X})
    return prediction

# Get probabilities for predictions
def predict_probability(X, parameters, axis=-1):
    # input -- X (dataset used to make prediction), papameters after training
    # output -- prediction
    W1 = tf.convert_to_tensor(parameters["W1"])
    b1 = tf.convert_to_tensor(parameters["b1"])
    W2 = tf.convert_to_tensor(parameters["W2"])
    b2 = tf.convert_to_tensor(parameters["b2"])
    W3 = tf.convert_to_tensor(parameters["W3"])
    b3 = tf.convert_to_tensor(parameters["b3"])
    W4 = tf.convert_to_tensor(parameters["W4"])
    b4 = tf.convert_to_tensor(parameters["b4"])
    params = {"W1": W1, "b1": b1, "W2": W2, "b2": b2, "W3": W3, "b3": b3, "W4": W4, "b4": b4}
    x = tf.placeholder("float")
    z4 = forward_propagation_for_predict(x, params)
    p = tf.nn.softmax(z4, axis=axis)
    sess = tf.Session()
    prediction = sess.run(p, feed_dict = {x: X})
    return prediction

# Build the model
def model(X_train, Y_train, starting_learning_rate = 0.0001, num_epochs = 1500, minibatch_size = 128, print_cost = True):
    # input -- X_train (training set), Y_train(training labels), X_test (test set), Y_test (test labels),
    # output -- trained parameters
    ops.reset_default_graph() # to be able to rerun the model without overwriting tf variables
    tf.set_random_seed(3)
    seed = 3
    (nf, ns) = X_train.shape
    nt = Y_train.shape[0]
    costs = []
    # create placeholders of shape (nf, nt)
    X, Y = create_placeholders(nf, nt)
    # initialize parameters
    parameters = initialize_parameters(nf=nf, ln1=100, ln2=50, ln3=25, nt=nt)
    # forward propagation: build the forward propagation in the tensorflow graph
    Z4 = forward_propagation(X, parameters)
    # cost function: add cost function to tensorflow graph
    cost = compute_cost(Z4, Y, parameters, 0.005)
    # Use learning rate decay
    global_step = tf.Variable(0, trainable=False)
    learning_rate = tf.train.exponential_decay(starting_learning_rate, global_step, 1000, 0.95, staircase=True)
    # backpropagation: define the tensorflow optimizer, AdamOptimizer is used.
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate)
    trainer = optimizer.minimize(cost, global_step=global_step)
    # initialize all the variables
    init = tf.global_variables_initializer()
    # start the session to compute the tensorflow graph
    with tf.Session() as sess:
        # run the initialization
        sess.run(init)
        # do the training loop
        for epoch in range(num_epochs):
            epoch_cost = 0.
            num_minibatches = int(ns / minibatch_size)
            seed = seed + 1
            minibatches = random_mini_batches(X_train, Y_train, minibatch_size, seed)
            for minibatch in minibatches:
                # select a minibatch
                (minibatch_X, minibatch_Y) = minibatch
                # run the session to execute the "optimizer" and the "cost", the feedict contains a minibatch for (X,Y).
                _ , minibatch_cost = sess.run([trainer, cost], feed_dict={X: minibatch_X, Y: minibatch_Y})
                epoch_cost += minibatch_cost / num_minibatches
            # print the cost every epoch
            if print_cost == True and (epoch+1) % 5 == 0:
                print ("Cost after epoch %i: %f" % (epoch+1, epoch_cost))
                costs.append(epoch_cost)
        parameters = sess.run(parameters)
        print ("Parameters have been trained!")
        # calculate the correct predictions
        correct_prediction = tf.equal(tf.argmax(Z4), tf.argmax(Y))
        # calculate accuracy on the test set
        accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
        print ("Train Accuracy:", accuracy.eval({X: X_train, Y: Y_train}))
        return parameters

# Function to train model [CLP 6/19/2020]
# - Inputs:
#    * train_set = pandas dataframe of expression data with genes as rows and cells as columns
#    * train_label = list of labels for each cell to train on
#    * learning_rate = the learning rate for the neural network
#    * num_epocs = the number of times to pass data through neural netowrk
#    * minibatch_size = size of batches the data is broken into
#    * print_cost = wether to print the cost for each five epochs
# - Outputs:
#    * parameters = the parameters from the trained neural network (the classifier)
#    * label_to_type_dict = a dictionary to convert the numbers from the classifier back into the original labels
#    * genes = the order and names of genes used to train the model
def train_model(train_set, train_label, learning_rate = 0.0001, num_epochs = 200, minibatch_size = 128, print_cost = True):
    nt = len(set(train_label))
    train_set, genes = scale_single_set(train_set, train = True)
    type_to_label_dict = dict(zip(set(train_label),range(len(set(train_label)))))
    label_to_type_dict = dict(zip(range(len(set(train_label))),set(train_label)))
    print("Cell Types in training set:", type_to_label_dict)
    print("# Training cells:", len(train_label))
    train_label = [type_to_label_dict[i] for i in train_label]
    train_label = one_hot_matrix(train_label, nt)
    parameters = model(train_set, train_label, learning_rate, num_epochs, minibatch_size, print_cost)
    return (parameters, label_to_type_dict, genes)

# Function to predict on new data [CLP 6/19/2020]
# - Inputs:
#    * new_data = pandas dataframe of expression data with genes as rows and cells as columns
#    * parameters = the parameters from the trained neural network (the classifier)
#    * label_to_type_dict = a dictionary to convert the numbers from the classifier back into the original labels
#    * genes = the order and names of genes used to train the model
# - Outputs:
#    * predicted_label = predicted labels for all cells
def predict_new_data(new_data, parameters, label_to_type_dict, genes):
    new_data = new_data.loc[genes]
    test_set = scale_single_set(new_data, train = False)
    test_predict = predict(test_set, parameters)
    predicted_label = pd.DataFrame({"cellname":new_data.columns, "celltype":[label_to_type_dict[i] for i in test_predict]})
    return predicted_label

# Function to predict probabilities [CLP 6/19/2020]
# - Inputs:
#    * new_data = pandas dataframe of expression data with genes as rows and cells as columns
#    * parameters = the parameters from the trained neural network (the classifier)
#    * label_to_type_dict = a dictionary to convert the numbers from the classifier back into the original labels
#    * genes = the order and names of genes used to train the model
# - Outputs:
#    * predicted_probabilities = probabilities for predicted labels for all cells
def predict_probabilities(new_data, parameters, genes, axis=-1):
    new_data = new_data.loc[genes]
    test_set = scale_single_set(new_data, train = False)
    test_probs = pd.DataFrame(predict_probability(test_set, parameters, axis=axis))
    #test_predict.index = [label_to_type_dict[x] for x in range(test_predict.shape[0])]
    #test_predict.columns = test_set.columns
    return test_probs

