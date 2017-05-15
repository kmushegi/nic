"""
Neural Networks for Digit Recognition - Project 4
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Project 5. This file 
contains the implementation of data loading routine used to read the training and
testing data on the cifar10 dataset.

*_data is a tuple. First element is array of arrays of size 1024 containing
the binary image; the second element is the array of arrays containing the
information about the "desired" output formatted in a 10 output neuron way.
[1,0,0,0,0,0,0,0,0,0] ==> 0
[0,0,1,0,0,0,0,0,0,0] ==> 2
[0,0,0,0,0,0,0,0,0,1] ==> 9 and etc.
"""

from __future__ import absolute_import
import cPickle
import keras
from keras.datasets import cifar10
import numpy as np
import sys

np.set_printoptions(threshold=np.nan)

bitmap_training_data_file = "/home/kmushegi/p5/data/bitmaps/optdigits-32x32.tra"
bitmap_testing_data_file = "/home/kmushegi/p5/data/bitmaps/optdigits-32x32.tes"
downsampled_training_data_file = "/home/kmushegi/p5/data/downsampled/optdigits-8x8-int.tra"
downsampled_testing_data_file = "/home/kmushegi/p5/data/downsampled/optdigits-8x8-int.tes"
cifar_dir = "/home/kmushegi/p5/data/cifar-10-batches-py/"

NUM_TRAINING_IMAGES = 3823

def read_data_down_sampled(data_file):
	f = open(data_file)

	#Containers for data and actual numbers
	inputs = []
	desired_outputs = []

	for line in f:
		tokens = line.split(',')

		#Each line in data file is an element
		inputs.append(" ".join(tokens[:-1]))
		#Actual number for corresponding data elements
		desired_outputs.append(int(tokens[-1]))

	return (inputs, desired_outputs)

def read_data_bit_map(data_file):
	f = open(data_file)

	for i in xrange(3): #skip first 3 lines
		next(f)

	#Containers bits and values
	inputs = []
	desired_outputs = []

	#Tracks each bitmap for values
	counter = 0

	current_input = ""
	for line in f:
		#While haven't reached end of bitmap
		if(counter < 32):
			current_input += line.strip()
			counter += 1
		#Insert this bitmap into input data
		else:
			inputs.append(current_input)
			desired_outputs.append(int(line))

			current_input = ""
			counter = 0

	return (inputs, desired_outputs)

def load_majercik_data(bit_map, num_output_neurons):
	if bit_map:
		#(inputs , desired_outputs) for training and testing data
		(tr_i, tr_o) = read_data_bit_map(bitmap_training_data_file)
		(te_i,te_o) = read_data_bit_map(bitmap_testing_data_file)

		#Change each 32x32 bitmap to a list of 1024 elements
		training_inputs = [np.reshape(list(x), (1024, 1)) for x in tr_i]
		testing_inputs = [np.reshape(list(x), (1024, 1)) for x in te_i]

		#Makes list of 10 elements with 1.0 at index = desired_output
		if num_output_neurons == 10:
			training_outputs = [digit_to_vector_representation(x) for x in tr_o]
			testing_outputs = [digit_to_vector_representation(x) for x in te_o]
		else:#Only 1 output node
			training_outputs = [int(x) for x in tr_o]
			testing_outputs = [int(x) for x in te_o]

		#Creates a list where each element is the 1024-list paired with corr. desired output
		training_data = zip(training_inputs,training_outputs)
		testing_data = zip(testing_inputs,testing_outputs)

		return (training_data, testing_data)
	else: #down-sampled
		#(inputs , desired_outputs) for training and testing data
		(tr_i, tr_o) = read_data_down_sampled(downsampled_training_data_file)
		(te_i, te_o) = read_data_down_sampled(downsampled_testing_data_file)

		#Change list of bits and spaces to a list of 64 bits
		training_inputs = [np.reshape(x.split(' '), (64, 1)) for x in tr_i]
		testing_inputs = [np.reshape(x.split(' '), (64, 1)) for x in te_i]

		#Makes list of 10 elements with 1.0 at index = desired_output
		if num_output_neurons == 10:
			training_outputs = [digit_to_vector_representation(x) for x in tr_o]
			testing_outputs = [digit_to_vector_representation(x) for x in te_o]
		else:#Only 1 output node
			training_outputs = [int(x) for x in tr_o]
			testing_outputs = [int(x) for x in te_o]

		#Creates a list where each element is the 64-bit list paired with corr. desired output
		training_data = zip(training_inputs, training_outputs)
		testing_data = zip(testing_inputs, testing_outputs)

		return (training_data, testing_data)

def load_data_batch(fp):
	f = open(fp,'rb')
	data_batch = cPickle.load(f)
	data_batch_dec = {}

	for key, value in data_batch.items():
		data_batch_dec[key.decode('utf8')] = value

	data_batch = data_batch_dec
	f.close()

	data = data_batch['data']
	labels = np.array(data_batch['labels'])

	return data,labels

def load_cifar_data_nn(n_batches):
	N_TRAINING_DATA = 50000

	for b in xrange(1,n_batches):
		batch_path = cifar_dir + "data_batch_" + str(b)
		print batch_path
		data, labels = load_data_batch(batch_path)

		if b == 1:
			training_inputs = data
			training_outputs = labels
		else:
			training_inputs = np.append(training_inputs,data,axis=0)
			training_outputs = np.append(training_outputs,labels,axis=0)

	training_inputs = normalize(training_inputs)
	training_inputs = training_inputs.reshape(-1,3072,1)
	training_outputs = [digit_to_vector_representation(y) for y in training_outputs]

	batch_path = cifar_dir + "test_batch"
	testing_inputs, testing_outputs = load_data_batch(batch_path)

	testing_inputs = normalize(testing_inputs)
	testing_inputs = testing_inputs.reshape(-1,3072,1)
	testing_outputs = [digit_to_vector_representation(y) for y in testing_outputs]

	training_data = zip(training_inputs, training_outputs)
	testing_data = zip(testing_inputs, testing_outputs)

	return (training_data, testing_data)

def load_cifar_data_cnn():
	(x_train,y_train),(x_test,y_test) = cifar10.load_data()
	y_train = keras.utils.to_categorical(y_train, 10)
	y_test = keras.utils.to_categorical(y_test, 10)
	x_train = x_train.astype('float32')
	x_test = x_test.astype('float32')
	x_train /= 255.0
	x_test /= 255.0

	return (x_train,y_train),(x_test,y_test)

def get_data(which_network,dataset='bitmap', num_output_neurons=10):
	if(dataset == "bitmap"):
		return load_majercik_data(bit_map=1,num_output_neurons=num_output_neurons)
	elif(dataset == "downsampled"):
		return load_majercik_data(bit_map=0,num_output_neurons=num_output_neurons)
	elif(dataset == "cifar10"):
		if(which_network == 'nn'):
			return load_cifar_data_nn(n_batches=3)
		elif(which_network == 'cnn'):
			return load_cifar_data_cnn()

#Creates 10 element list of 0.0's except with 1.0 at jth index
def digit_to_vector_representation(j):
	e = np.zeros((10,1))
	e[j] = 1.0
	return e

def normalize(data):
	return data/255.0










