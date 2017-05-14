"""
Neural Networks for Image Classification - Final Project
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Neural Networks for Image Classification, Final Project. This file
contains the implementation of NEED TO WRITE THIS
"""

import sys
import random
import math
import network as nn
import cnetwork as cnn

class GA(object):

	def __init__(self, networkType, numIterations, populationSize, selection, crossover, crossoverProb, mutationProb, nnParams):
		self.networkType = networkType
		self.numberOfIterations = numIterations
		self.populationSize = populationSize
		self.selection = positiveLearningRate
		self.crossover = crossover
		self.crossoverProb = crossoverProb
		self.mutationProb = mutationProb
		self.nnParams = nnParams
		self.population = initializePopulation()

	def initializePopulation(self):
		pop = []
		for _ in range(0, self.populationSize):
			if (self.networkType == "nn"):
				temp = []
				for key in nnParams:
					temp.append(random.choice(self.nnParams[key]))
				net = nn.Network(temp[0] , temp[1] , temp[2] , temp[3] ) #Parameters
			else:
				temp = []
				for key in nnParams:
					temp.append(random.choice(self.nnParams[key]))
				net = nn.Network(temp[0] , temp[1] , temp[2] , temp[3] , temp[4] , temp[5] , temp[6] , temp[7] ,temp[8] , temp[9]) #Parameters
		pop.append(net)
		return pop

	def evalFitness(network):



	def runGA(self):


