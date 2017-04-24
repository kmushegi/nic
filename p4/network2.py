import numpy as np

class Network(object):

	#sizes = array with number of neurons in each level: input / hidden / output
	#[32,15,10] -> 32 input, 15 hidden and 10 output neurons
	def __init__(self, sizes):
		self.num_layers = len(size)
		self.sizes = sizes
		self.biases = [np.random.randn(y,1) for y in sizes[1:]]
		self.weights = [np.random.randn(y,x) for x, y in zip(sizes[:-1],sizes[1:])]


def sigmoid(z):
	return 1.0/(1.0+np.exp(-z))

def sigmoidDerivative(z):
	return sigmoid(z)*(1-sigmoid(z))