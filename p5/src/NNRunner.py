"""
Neural Networks for Digit Recognition - Project 4
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Neural Networks for Digit Recognition, Project 4. This file
is an entry point into the program. It is able to take parameters for NN parameters
from the command line if CL_PARAMS is set to True. This file also reads the 
training and testing data and spawns the neural network training process.
"""

import data_loader as dl
import network
import sys

N_EPOCHS = 50

#if # of CL parameters is correct take them, else go to default
if len(sys.argv) == 4:
	bit_map = int(sys.argv[1]) #1 - bitmap, 0 - downsampled
	n_hid_neurons = int(sys.argv[2]) #hidden layer nodes
	n_out_neurons = int(sys.argv[3]) # 10 or 1
	l_rate = float(sys.argv[4])
else:
	bit_map = 1 # 1 - bitmap, 0 - downsampled
	n_hid_neurons = 30
	n_out_neurons = 10 # 10 or 1
	l_rate = 0.5

assert(n_out_neurons == 10 or n_out_neurons == 1)

if bit_map:
	n_in_neurons = 1024
else:
	n_in_neurons = 64

#obtain data using data_loader.py
(training_data, test_data) = dl.get_data("bitmap",n_out_neurons)

#create an instance of the neural network with specified parameters
net = network.Network(num_inputs=n_in_neurons, num_hidden=n_hid_neurons, num_outputs=n_out_neurons, 
						num_epochs=N_EPOCHS, learning_rate=l_rate)

#start the training process
net.train(training_data, test_data)

