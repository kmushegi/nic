"""
Neural Networks for Digit Recognition - Project 4
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part Project 4. This file is an entry point into the Perceptron program. 
It is able to take parameters for the NN if correct # of them is specified. This file also reads the 
training and testing data and spawns the neural network training process.
"""

import data_loader as dl
import network
import sys

N_EPOCHS = 50

#if # of CL parameters is correct take them, else go to default
if len(sys.argv) == 4:
	dataset = sys.argv[1] #bitmap, downsampled, cifar10
	n_out_neurons = int(sys.argv[3]) # 10 or 1
	l_rate = float(sys.argv[4])
else:
	dataset = "bitmap" #bitmap,downsampled, cifar10
	n_hid_neurons = 100
	n_out_neurons = 10 # 10 or 1
	l_rate = 0.5

assert(n_out_neurons == 10 or n_out_neurons == 1)

if(dataset == "bitmap"):
	n_in_neurons = 1024
elif(dataset == "downsampled"):
	n_in_neurons = 64
elif(dataset == "cifar10"):
	n_in_neurons = 3072
else:
	print("Choose dataset [bitmap, downsampled, cifar10]")
	sys.exit(1)

#obtain data using data_loader.py
(training_data, test_data) = dl.get_data(dataset,n_out_neurons)

#create an instance of the neural network with specified parameters
net = network.Network([n_in_neurons,100,50,n_out_neurons],n_epochs=N_EPOCHS, lr=l_rate)

#start the training process
net.train(training_data, test_data)

