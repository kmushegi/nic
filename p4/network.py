import numpy as np
import sys
import math

#np.set_printoptions(threshold=np.nan)

class Network(object):
	#sizes = array with number of neurons in each level: input / hidden / output
	#[32,15,10] -> 32 input, 15 hidden and 10 output neurons
	def __init__(self, numInputNodes, numOutputNodes, numEpochs, learningRate):
		self.numInputNodes = numInputNodes + 1 #Add 1 for bias node
		self.numOutputNodes = numOutputNodes
		self.numEpochs = numEpochs
		self.learningRate = learningRate
		self.inputNodes = [1.0] * self.numInputNodes
		self.outputNodes = [1.0] * self.numOutputNodes
		self.weights = np.random.rand(self.numInputNodes, self.numOutputNodes) * 0.15

	def train(self, training_data):
		for _ in range(self.numEpochs):
			for trainingSet in training_data:
				self.initializeInputNodes(trainingSet[0])
				self.calculateOutputNodeValues()
				self.updateWeights(trainingSet[1])

	def initializeInputNodes(self, trainingInput):
		for i in range(trainingInput.size):
			self.inputNodes[i] = 1.0 * int(trainingInput.item(i))
		self.inputNodes[-1] = 1.0 #Bias node

	def calculateOutputNodeValues(self):
		for i in range(len(self.outputNodes)):
			self.outputNodes[i] = self.sigmoidMod(self.sumInputs(i))
			print(self.outputNodes[i])

	def sumInputs(self, index):
		inputSum = 0.0
		for i in range(len(self.inputNodes)):
			inputSum += self.inputNodes[i] * self.weights[i][index]
		return inputSum

	def sigmoidMod(self, x):
		return 1.0/(1.0+math.exp(-x+0.5))

	def sigmoidDerivativeMod(self, x):
		return self.sigmoidMod(x)*(1-self.sigmoidMod(x))

	def sigmoid(self, x):
		return 1.0/(1.0+np.exp(-x))

	def updateWeights(self, expectedOutput):
		for i in range(self.numInputNodes):
			for j in range(self.numOutputNodes):
				self.weights[i][j] += self.weightUpdate(self.inputNodes[i], self.calculateError(expectedOutput, self.outputNodes), self.sumInputs(j))

	def weightUpdate(self, activationLevel, error, inputSum):
		return self.learningRate * activationLevel * error * self.sigmoidDerivativeMod(inputSum)

	def sigmoidDerivative(self, x):
		return self.sigmoid(x)*(1-self.sigmoid(x))

	def calculateError(self, expectedResult, outputResult):
			difference = [x1 - x2 for (x1, x2) in zip(expectedResult, outputResult)]
			return np.linalg.norm(difference)




