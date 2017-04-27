from __future__ import print_function
import numpy as np
import sys
import math
import random
import time

np.set_printoptions(threshold=10)

class Network(object):
	#sizes = array with number of neurons in each level: input / hidden / output
	#[32,15,10] -> 32 input, 15 hidden and 10 output neurons
	def __init__(self, num_inputs, num_outputs, num_epochs, learning_rate):
		self.num_inputs = num_inputs + 1 #Add 1 for bias node
		self.num_outputs = num_outputs
		self.num_epochs = num_epochs
		self.learning_rate = learning_rate

		self.in_activations = np.reshape(np.ones(self.num_inputs),(self.num_inputs,1))
		self.out_activations = np.ones(self.num_outputs) #[1.0] * self.num_outputs

		#Negative weights between -0.15 and 0.15
		self.weights = 0.3 * np.random.rand(self.num_inputs, self.num_outputs) - 0.15
		self.inputSums = np.ones(self.num_outputs) #[1.0] * self.num_outputs #REVIEW

		self.delta_io = np.zeros((self.num_inputs, self.num_outputs))

	def sigmoid(self, x):
		return 1.0 / (1.0 + np.exp(-x))

	def sigmoidPrime(self, x):
		return self.sigmoid(x) * (1.0 - self.sigmoid(x))

	def train(self, training_data, test_data=None):

		if test_data: n_test_samples = len(test_data)

		for i in xrange(self.num_epochs):
			random.shuffle(training_data)
			for index, sample in enumerate(training_data):
				print("Sample: ",index,end='\r')
				sys.stdout.flush()
				self.feedForward(sample[0])
				error = self.update(sample[1])
			print("")
			if test_data:
				print("Epoch {0}: {1} / {2} = {3}%".format(
					i, self.test(test_data),n_test_samples))
			else:
				print("Epoch {0}".format(i))

	def feedForward(self, inputs):
		self.in_activations[0:self.num_inputs-1] = inputs
		self.inputSums = np.dot(self.weights.T, self.in_activations)
		self.out_activations = self.sigmoid(self.inputSums)
		return self.out_activations

	def update(self, desired_output):
		error = -(desired_output - self.out_activations)
		out_deltas = self.sigmoidPrime(self.out_activations) * error

		#update weights from input to output layer
		for i in xrange(self.num_inputs):
			for o in xrange(self.num_outputs):
				delta = out_deltas[o] * self.in_activations[i]
				self.weights[i][o] -= self.learning_rate * delta + self.delta_io[i][o]
				self.delta_io[i][o] = delta

		error = 0.0
		for do in xrange(len(desired_output)):
			error += 0.5 * (desired_output[do] - self.out_activations[do]) ** 2
		return error
		'''
		for i in xrange(self.num_inputs):
			for j in xrange(self.num_outputs):
				self.weights[i][j] += self.learning_rate * self.in_activations[i] * (expectedOutput[j] - self.out_activations[j]) * self.sigmoidPrime(self.inputSums[j])
		'''

	def test(self, test_data):

		test_results = [(np.argmax(self.feedForward(x)), np.argmax(y))
						for (x,y) in test_data]
		return sum((x == y) for (x,y) in test_results)


