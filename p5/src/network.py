"""
Evolving Neural Networks Using Genetic Algorithms - Project 5
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Project 5. This file contains the implementation of the Multi-layer
Perceptron model of artificial neural networks with an arbitrary number of hidden layers.
"""

from __future__ import print_function
import numpy as np
import sys
import random
import math

np.set_printoptions(threshold=np.nan)
np.seterr(all='ignore')


#momentum term - magic number selected by trial and error, let's us go further
#towards gradient descent, but not too far to miss the minimum point
M = 0.5

class Network(object):

	def __init__(self, layer_info, n_epochs, lr, startWeightRange): #layer_info => [in,hidden,...,hidden,out] 
		self.n_layers = len(layer_info)
		self.layer_info = layer_info
		self.n_inputs = layer_info[0] + 1 #add 1 for bias node
		self.n_epochs = n_epochs
		self.lr = lr
		self.parameters = [layer_info, n_epochs, lr, startWeightRange]

		self.weights = []		#weights[i] = weights between layer i and layer i+1
		self.deltas = []		#deltas[i] = deltas between layer i and layer i+1
		self.activations = []	#activations[i] = activations of layer i

		#unlike everything else there are (n_layers) activations
		in_activations = np.reshape(np.ones(self.n_inputs),(self.n_inputs,1))
		self.activations.append(in_activations)

		#initialize random weights in [-0.15,0.15], deltas for the momentum term,layer activations
		for i in xrange(self.n_layers - 1):
			if i == 0:
				temp_w = (2*startWeightRange) * np.random.randn(self.n_inputs, self.layer_info[i+1]) - startWeightRange
				temp_d = np.zeros((self.n_inputs, self.layer_info[i+1]))
			else:
				temp_w = (2*startWeightRange) * np.random.randn(self.layer_info[i], self.layer_info[i+1]) - startWeightRange
				temp_d = np.zeros((self.layer_info[i], self.layer_info[i+1]))

			#append generated weights, momentum deltas and layer activations to global containers
			self.weights.append(temp_w)
			self.deltas.append(temp_d)
			self.activations.append(np.ones(self.layer_info[i]))

	#activation function
	def sigmoid(self, x):
		return 1.0 / (1.0 + np.exp(np.longfloat(-x)))

	#derivative of activation function
	def sigmoidPrime(self, x):
		return self.sigmoid(x) * (1.0 - self.sigmoid(x))

	#train the neural network on training_data with optional test_data for evaluation
	def train(self, training_data, test_data=None):
		if test_data: n_test_samples = len(test_data)

		for i in xrange(self.n_epochs):
			random.shuffle(training_data) #shuffle training data for 'random' sampling
			
			for index, sample in enumerate(training_data): #for every input/output pair
				print("Sample: ",index,end='\r')
				sys.stdout.flush()
				self.feedForward(sample[0])
				self.update(sample[1])
			print("")
			if test_data:
				test_results = self.test(test_data)
				print("Epoch {0}:\t {1} \t {2}".format(i+1, test_results,n_test_samples))
			else:
				print("Epoch {0}".format(i+1))

		return test_results

	#run the input through the neural network with existing weights
	def feedForward(self, inputs):
		#set input layer nodes to actual input, except for the bias node
		self.activations[0][0:self.n_inputs-1] = inputs

		#feed the input through the neural network input ---> output
		for i in xrange(self.n_layers - 1):
			sum = np.dot(self.weights[i].T,self.activations[i])
			self.activations[i+1] = self.sigmoid(sum)

		return self.activations[-1] #return activation of the output layer

	#update weights in the direction of the gradient descent
	def update(self, desired_output):

		t_delta = None

		#update the weights starting at the output layer and moving backwards -- 'backpropagate'
		for i in xrange(1,self.n_layers):
			if i == 1: #error term = desired - out_activation at the last layer
				error = -(desired_output - self.activations[-i])
				t_delta = self.sigmoidPrime(self.activations[-i]) * error
			else:
				error = np.dot(self.weights[-(i-1)],t_delta) #compute error
				t_delta = self.sigmoidPrime(self.activations[-i]) * error #compute t_delta

			delta = t_delta.T * self.activations[-(i+1)] #compute the direction of weight change
			self.weights[-i] -= (self.lr * delta + M * self.deltas[-i]) #update weights
			self.deltas[-i] = delta #store the current delta for next epoch momentum term

	#feed test_data through the neural net and return the number of correct predictions
	def test(self, test_data):
		if self.layer_info[-1] == 10:
			#'zip' together perceptron results and the expected test output to form a tuple
			test_results = [(np.argmax(self.feedForward(x)), np.argmax(y)) for (x,y) in test_data]
		elif self.layer_info[-1] == 1:
			#perceptron output is multiplied by 10, as our expected test output is a number [0,9]
			test_results = [(math.floor(self.feedForward(x) * 10.0), y) for (x,y) in test_data]
		#return the total number of test cases that were correctly classified
		return sum((x == y) for (x,y) in test_results)
		








