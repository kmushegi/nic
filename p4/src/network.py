"""
Neural Networks for Digit Recognition - Project 4
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Neural Networks for Digit Recognition, Project 4. This file
contains the implementation of the Perceptron model of artificial neural 
networks, with two layers: input & output.
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

	def __init__(self, num_inputs, num_outputs, num_epochs, learning_rate):
		self.num_inputs = num_inputs + 1 #add 1 for bias node
		self.num_outputs = num_outputs
		self.num_epochs = num_epochs
		self.learning_rate = learning_rate

		self.in_activations = np.reshape(np.ones(self.num_inputs),(self.num_inputs,1))
		self.out_activations = np.ones(self.num_outputs)

		#keep track of how much weights need to change at next iteration
		self.delta_io = np.zeros((self.num_inputs, self.num_outputs))

		#initialize weights randomly between -0.15 and 0.15
		self.weights = 0.3 * np.random.rand(self.num_inputs, self.num_outputs) - 0.15
		self.inputSums = np.ones(self.num_outputs)

	#activation function
	def sigmoid(self, x):
		return 1.0 / (1.0 + np.exp(np.longfloat(-x)))

	#derivative of activation function
	def sigmoidPrime(self, x):
		return self.sigmoid(x) * (1.0 - self.sigmoid(x))

	#start the training process, test_data is optional and is used for evaluation
	def train(self, training_data, test_data=None):
		if test_data: n_test_samples = len(test_data)

		#for every epoch
		for i in xrange(self.num_epochs):
			random.shuffle(training_data) #shuffle training data for 'random' sampling
			
			for index, sample in enumerate(training_data): #for every input/output pair
				self.feedForward(sample[0])
				self.update(sample[1])
			if test_data:
				print("Epoch {0}:\t {1} \t {2}".format(i+1, self.test(test_data),n_test_samples))
			else:
				print("Epoch {0}".format(i+1))

		return self.weights

	#run the input through the neural network with existing weights
	def feedForward(self, inputs):
		#set input nodes to actual input, except for the bias node
		self.in_activations[0:self.num_inputs-1] = inputs

		#dot weights with the input
		self.inputSums = np.dot(self.weights.T, self.in_activations)

		#Plug the sum into the activtion function
		self.out_activations = self.sigmoid(self.inputSums)
		return self.out_activations

	#update weights in the direction of the gradient descent
	def update(self, desired_output):

		if (self.num_outputs == 1): #'normalize' desired output to be in [0,1]
			desired_output /= 10.0

		#calculate error between current perceptron output and expected output
		error = -(desired_output - self.out_activations)
		out_deltas = self.sigmoidPrime(self.out_activations) * error

		delta = out_deltas.T * self.in_activations
		self.weights -= (self.learning_rate * delta + M * self.delta_io)
		self.delta_io = delta

	#feed test_data through the neural net and return the number of correct predictions
	def test(self, test_data):
		if self.num_outputs == 10:
			#'zip' together perceptron results and the expected test output to form a tuple
			test_results = [(np.argmax(self.feedForward(x)), np.argmax(y)) for (x,y) in test_data]
		elif self.num_outputs == 1:
			#perceptron output is multiplied by 10, as our expected test output is a number [0,9]
			test_results = [(math.floor(self.feedForward(x) * 10.0), y) for (x,y) in test_data]
		#return the total number of test cases that were correctly classified
		return sum((x == y) for (x,y) in test_results)







