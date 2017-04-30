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

CL_PARAMS = False

if CL_PARAMS:
	bit_map = int(sys.argv[1]) #1 - bitmap, 0 - downsampled
	n_out_neurons = int(sys.argv[2]) # 10 or 1
	n_epochs = int(sys.argv[3])
	l_rate = float(sys.argv[4])
else:
	bit_map = 1 # 1 - bitmap, 0 - downsampled
	n_out_neurons = 1 # 10 or 1
	n_epochs = 10
	l_rate = 0.01

assert(n_out_neurons == 10 or n_out_neurons == 1)

if bit_map:
	n_in_neurons = 1024
else:
	n_in_neurons = 64


(training_data, test_data) = dl.format_data(bit_map,n_out_neurons)

net = network.Network(num_inputs=n_in_neurons, num_outputs=n_out_neurons, 
						num_epochs=n_epochs, learning_rate=l_rate)

net.train(training_data, test_data)

