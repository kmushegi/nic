from __future__ import print_function
import numpy as np
import sys
import math
import random
import time

np.set_printoptions(threshold=np.nan)

class Network(object):
	#sizes = array with number of neurons in each level: input / hidden / output
	#[32,15,10] -> 32 input, 15 hidden and 10 output neurons
	def __init__(self, num_inputs, num_outputs, num_epochs, learning_rate):
		self.num_inputs = num_inputs + 1 #Add 1 for bias node
		self.num_outputs = num_outputs
		self.num_epochs = num_epochs
		self.learning_rate = learning_rate

		self.in_activations = np.reshape(np.ones(self.num_inputs),(self.num_inputs,1))
		self.out_activations = [1.0] * self.num_outputs

		#Negative weights between -0.15 and 0.15
		self.weights = 0.3 * np.random.rand(self.num_inputs, self.num_outputs) - 0.15
		self.inputSums = [1.0] * self.num_outputs

	def train(self, training_data, test_data=None):
		for i in xrange(self.num_epochs):
			random.shuffle(training_data)
			for index, sample in enumerate(training_data):
				print("Sample: " + str(index), end='\r')
				sys.stdout.flush()
				stf = time.time()

				self.feedForward(sample[0])

				ef = time.time()

				self.update(sample[1])

				ew = time.time()

				print("Update Weights: " + str(ew-ef) + "\t FF: " + str(ef-stf))
			# print "Epoch: {0} done.".format(i)
		if test_data:
			self.evaluate(test_data)

	def evaluate(self, test_data):
		results = [(np.argmax(self.feedForward(x)), y) for x,y in test_data]
		#print(results)

	def feedForward(self, inputs):
		self.in_activations[0:self.num_inputs-1] = inputs
		self.inputSums = np.dot(self.weights.T, self.in_activations)
		self.out_activations = self.sigmoid(self.inputSums)

	def sigmoid(self, x):
		return 1.0 / (1.0 + np.exp(-x))

	def sigmoidPrime(self, x):
		return self.sigmoid(x) * (1.0 - self.sigmoid(x))

	def update(self, desired_output):
		error = -(desired_output - self.out_activations)
		out_deltas = self.sigmoidPrime(self.out_activations) * error

		#update weights from input to output layer
		for i in xrange(self.num_inputs):
			for o in xrange(self.num_outputs):
				delta = out_deltas[o] * self.in_activations[i]
				self.weights[i][o] -= 





		for i in xrange(self.num_inputs):
			for j in xrange(self.num_outputs):
				self.weights[i][j] += self.learning_rate * self.in_activations[i] * (expectedOutput[j] - self.out_activations[j]) * self.sigmoidPrime(self.inputSums[j])




