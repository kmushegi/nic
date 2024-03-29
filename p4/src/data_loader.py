"""
Neural Networks for Digit Recognition - Project 4
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Neural Networks for Digit Recognition, Project 4. This file 
contains the implementation of data loading routine used to read the training and
testing data from .tes and .tra files provided by Professor Majercik.
"""

import numpy as np

bitmap_training_data_file = "../data/bitmaps/optdigits-32x32.tra"
bitmap_testing_data_file = "../data/bitmaps/optdigits-32x32.tes"

downsampled_training_data_file = "../data/downsampled/optdigits-8x8-int.tra"
downsampled_testing_data_file = "../data/downsampled/optdigits-8x8-int.tes"

NUM_TRAINING_IMAGES = 3823

np.set_printoptions(threshold=np.nan)

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

'''
*_data is a tuple. First element is array of arrays of size 1024 containing
the binary image; the second element is the array of arrays containing the
information about the "desired" output formatted in a 10 output neuron way.
[1,0,0,0,0,0,0,0,0,0] ==> 0
[0,0,1,0,0,0,0,0,0,0] ==> 2
[0,0,0,0,0,0,0,0,0,1] ==> 9 and etc.
'''

def format_data(bit_map, num_output_neurons):
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

#Creares 10 element list of 0.0's except with 1.0 at jth index
def digit_to_vector_representation(j):
	e = np.zeros((10,1))
	e[j] = 1.0
	return e


#(tr_data, test_data) = format_data(0,10)
#print(np.array_str(tr_data[0][0]))
#assert(len(tr_data) == NUM_TRAINING_IMAGES)
