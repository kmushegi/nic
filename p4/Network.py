import numpy as np
import sys

#np.set_printoptions(threshold=np.nan)

class Network(object):

	#sizes = array with number of neurons in each level: input / hidden / output
	#[32,15,10] -> 32 input, 15 hidden and 10 output neurons
	def __init__(self, numInputNodes, numOutputNodes, numEpochs, learningRate):
		self.numInputNodes = numInputNodes + 1 #Add 1 for bias node
		self.numOutputNodes = numOutputNodes
		self.numEpochs = numEpochs
		self.learningRate = learningRate
		self.inputNodes = [1.0] * numInputNodes
		self.outputNodes = [1.0] * numOutputNodes
		self.weights = 0.15 * np.random.rand(self.numInputNodes, self.numOutputNodes)

	def train(self, training_data):

		# sys.exit(1)

		for _ in range(self.numEpochs):
			for trainingSet in training_data:
				self.initializeInputNodes(trainingSet[0])
				self.calculateOutputNodeValues()
				self.updateWeights()
				sys.exit(1)

	def updateWeights(self):
		return 1

	def calculateOutputNodeValues(self):
		for i in range(len(self.outputNodes)):
			self.outputNodes[i] = self.sigmoid(self.sumInputs(i))

	def sumInputs(self, index):
		inputSum = 0.0
		for i in range(len(self.inputNodes)):
			print(self.inputNodes[i] * self.weights[i][index])
			inputSum += self.inputNodes[i] * self.weights[i][index]

		return inputSum

	def initializeInputNodes(self, trainingInput):
		for i in range(trainingInput.size):
			self.inputNodes[i] = 1.0 * int(trainingInput.item(i))

	def sigmoidMod(self, x):
		return 1.0/(1.0+np.exp(-x+0.5))

	def sigmoidDerivativeMod(x):
		return sigmoidMod(x)*(1-sigmoidMod(x))

	def sigmoid(self, x):
		return 1.0/(1.0+np.exp(-x))

	def sigmoidDerivative(x):
		return sigmoid(x)*(1-sigmoid(x))

	def calculateError(expectedResult, outputResult):
			difference = [x1 - x2 for (x1, x2) in zip(expectedResult, outputResult)]
			return np.sqrt(difference.dot(difference))

	def weightUpdate(learningRate, activationLevel, error, inputSum):
		return learningRate * activationLevel * error * sigmoidDerivative(inputSum)




