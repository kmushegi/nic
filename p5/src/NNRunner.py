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

which_network = 'nn' #or 'cnn'
nn_dataset = "cifar10" #bitmap,downsampled, cifar10
N_EPOCHS = 50

if which_network == 'nn':

	#if # of CL parameters is correct take them, else go to default
	if len(sys.argv) >= 5:
		dataset = sys.argv[1] #bitmap, downsampled, cifar10
		n_out_neurons = int(sys.argv[2]) # 10 or 1
		l_rate = float(sys.argv[3])
		startWeightR = float(sys.argv[4])

		hidden_layer_info = []
		for i in xrange(5,len(sys.argv)):
			hidden_layer_info.append(int(sys.argv[i]))
	else:
		dataset = 'cifar10' #bitmap,downsampled, cifar10
		n_out_neurons = 10 # 10 or 1
		l_rate = 1
		hidden_layer_info = [100,50]
		startWeightR = 0.15

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

	#load training & test data using data_loader.py
	(train_data, test_data) = dl.get_data(which_network,dataset,n_out_neurons)

	layer_info = [n_in_neurons]
	layer_info.extend(hidden_layer_info)
	layer_info.append(n_out_neurons)
	
	#create an instance of the neural network with specified parameters
	net = nn.Network(layer_info=layer_info, n_epochs=N_EPOCHS, lr=l_rate, startWeightRange=startWeightR)

	#start the training process
	net.train(train_data, test_data)

elif which_network == 'cnn':

	#if # of CL parameters is correct take them, else go to default
	if len(sys.argv) == 8:
		n_layers = int(sys.argv[1])
		dropout = int(sys.argv[2])
		batch_size = int(sys.argv[3])
		optimizer = sys.argv[4]
		data_augmentation = int(sys.argv[5])
		convActivation = sys.argv[6]
		denseActivation = sys.argv[7]
	else:
		n_layers = 4
		dropout = True
		batch_size = 32
		optimizer = 'sgd'
		data_augmentation = True
		convActivation = 'relu'
		denseActivation = 'softmax'

	(x_train,y_train),(x_test,y_test) = dl.get_data(which_network,'cifar10')

	net = cnn.CNetwork(n_epochs=N_EPOCHS,n_layers=n_layers,dropout=dropout,batch_size=batch_size,
		optimizer=optimizer,data_augmentation=data_augmentation,convActivation=convActivation,
		denseActivation=denseActivation,x_train=x_train,y_train=y_train,x_test=x_test,y_test=y_test)

	net.build_network()
	net.train()





















