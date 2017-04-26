import numpy as np
import sys
import math

#np.set_printoptions(threshold=np.nan)

class Network(object):
	#sizes = array with number of neurons in each level: input / hidden / output
	#[32,15,10] -> 32 input, 15 hidden and 10 output neurons
	def __init__(self, num_inputs, num_outputs, num_epochs, learning_rate):
		self.num_inputs = num_inputs + 1 #Add 1 for bias node
		self.num_outputs = num_outputs
		self.num_epochs = num_epochs
		self.learning_rate = learning_rate
		self.in_activations = [1.0] * self.num_inputs
		self.out_activations = [1.0] * self.num_outputs
		#Negative weights between -0.15 and 0.15
		self.weights = 0.3 * np.random.rand(self.num_inputs, self.num_outputs) - 0.15
		self.inputSums = [1.0] * self.num_outputs

	def train(self, training_data, test_data=None):
		for i in xrange(self.num_epochs):
			j = 1
			for trainingSet in training_data:
				print("New Training " + str(j))
				self.initialize_in_activations(trainingSet[0])
				self.calculateOutputNodeValues()
				self.updateWeights(trainingSet[1])
				j+=1
			print "Epoch: {0} done.".format(i)
		if test_data:
			self.evaluate(test_data)

	def evaluate(self, test_data):
		results = [(np.argmax(self.feedForward(x)), y) for x,y in test_data]
		#print(results)

	def feedForward(self, inputs):
		output = self.sigmoid(self.sumInputs(i))
		return output

	def initialize_in_activations(self, trainingInput):
		for i in range(trainingInput.size):
			self.in_activations[i] = 1.0 * int(trainingInput.item(i))
		self.in_activations[-1] = 1.0 #Bias node

	def calculateOutputNodeValues(self):
		for i in range(len(self.out_activations)):
			self.inputSums[i] = self.sumInputs(i)
			self.out_activations[i] = self.sigmoid(self.inputSums[i])
			print(self.out_activations[i])

	def sumInputs(self, index):
		inputSum = 0.0
		for i in range(len(self.in_activations)):
			inputSum += self.in_activations[i] * self.weights[i][index]
		return inputSum

	def sigmoid(self, x):
		return 1.0 / (1.0 + np.exp(-x))

	def sigmoidPrime(self, x):
		return self.sigmoid(x) * (1.0 - self.sigmoid(x))

	def updateWeights(self, expectedOutput):
		for i in range(self.num_inputs):
			for j in range(self.num_outputs):
				self.weights[i][j] += self.weightUpdate(self.in_activations[i], 
				expectedOutput[j] - self.out_activations[j], self.inputSums[j])

	def weightUpdate(self, activationLevel, error, inputSum):
		return self.learning_rate * activationLevel * error * self.sigmoidPrime(inputSum)




