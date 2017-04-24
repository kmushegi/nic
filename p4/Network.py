import numpy as np
import sys

class Network(object):

	#sizes = array with number of neurons in each level: input / hidden / output
	#[32,15,10] -> 32 input, 15 hidden and 10 output neurons
	def __init__(self, sizes):
		self.num_layers = len(sizes)
		self.sizes = sizes
		self.biases = [np.random.randn(y,1) for y in sizes[1:]]
		print("Biases size: " + str(len(self.biases)))
		self.weights = [np.random.randn(y,x) for x, y in zip(sizes[:-1],sizes[1:])]
		print("Weights size: " + str(len(self.weights)))

	def feedForward(self, i):
		for b, w in zip(self.biases, self.weights):
			a = sigmoid(np.dot(w, a) + b)
		return a

	def squaredGradientDescent(self, training_data, epochs, sample_size, learningRate):
		n = len(training_data)

		for e in xrange(epochs):
			for sample in training_data:
				self.process_sample(sample, learningRate)

	def process_sample(self, sample, learningRate):
		sums = self.biases
		#sums += (sample * self.weights)

def sigmoid(z):
	return 1.0/(1.0+np.exp(-z))

def sigmoidDerivative(z):
	return sigmoid(z)*(1-sigmoid(z))


