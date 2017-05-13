"""
Neural Networks for Digit Recognition - Project 5
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian
"""

import data_loader as dl
import network as nn
import cnetwork as cnn
import sys

which_network = 'cnn' #of 'cnn'
nn_dataset = "cifar10" #bitmap,downsampled, cifar10
N_EPOCHS = 50

if which_network == 'nn':

	#if # of CL parameters is correct take them, else go to default
	if len(sys.argv) == 4:
		dataset = sys.argv[1] #bitmap, downsampled, cifar10
		n_hid_neurons = int(sys.argv[2]) #hidden layer nodes
		n_out_neurons = int(sys.argv[3]) # 10 or 1
		l_rate = float(sys.argv[4])
	else:
		dataset = "cifar10" #bitmap,downsampled, cifar10
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

	#load training & test data using data_loader.py
	(train_data, test_data) = dl.get_data(dataset,n_out_neurons)

	#create an instance of the neural network with specified parameters
	net = nn.Network(num_inputs=n_in_neurons, num_hidden=n_hid_neurons, 
		num_outputs=n_out_neurons, num_epochs=N_EPOCHS, learning_rate=l_rate)

	#start the training process
	net.train(training_data, test_data)

elif which_network == 'cnn':
	(x_train,y_train),(x_test,y_test) = dl.get_data(which_network,'cifar10')
	n_classes = 10

	net = cnn.CNetwork(n_epochs=N_EPOCHS,n_layers=4,dropout=True,batch_size=32,optimizer='sgd',
		data_augmentation=True,convActivation='relu',denseActivation='softmax',
		x_train=x_train,y_train=y_train,x_test=x_test,y_test=y_test)

	net.build_network()
	net.train()





















