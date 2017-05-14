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
<<<<<<< Updated upstream
import rankedIndividual
=======
import data_loader as dl
>>>>>>> Stashed changes

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

		if (self.networkType == "nn"):
			self.train_dataANDtest_data = dl.get_data(which_network,dataset,n_out_neurons)

		else:
			self.x_train,y_trainANDx_test,y_test = dl.get_data(which_network,'cifar10')

	def initializePopulation(self):
		pop = []
		for _ in range(0, self.populationSize):
			if (self.networkType == "nn"):
<<<<<<< Updated upstream
				temp = []
				for key in nnParams:
					temp.append(random.choice(self.nnParams[key]))
				net = nn.Network(temp[0] , temp[1] , temp[2] , temp[3] , temp) #Parameters
=======
				n_in_neurons = 1024
				n_out_neurons = 10
				hidden_layer_info = random.choice(self.nnParams["hiddenInfo"])

				layer_info = [n_in_neurons]
				layer_info.extend(hidden_layer_info)
				layer_info.append(n_out_neurons)

				net = nn.Network(layer_info, random.choice(self.nnParams["epochs"]), random.choice(self.nnParams["learningRate"])) #Parameters

>>>>>>> Stashed changes
			else:
				temp = []
				for key in nnParams:
					temp.append(random.choice(self.nnParams[key]))
<<<<<<< Updated upstream
				net = cnn.CNetwork(temp[0] , temp[1] , temp[2] , temp[3] , temp[4] , temp[5] , temp[6] , temp[7] ,temp[8] , temp[9] , temp) #Parameters
		pop.append(net)
		return pop

	#def evalFitness(network):



	#def runGA(self):

=======
				net = nn.Network(temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7],) #Parameters

		pop.append(net)
		return pop

	# def runGA(self):
>>>>>>> Stashed changes

	def rs(self , pop):
		#Willbe the selected networks
		selected = []

		#Willbe used to calculate probabilites
		totalRank = ((len(selected)*(len(selected)+1))/2)

		#Will be a list of ranked individuals
		fitnessEvaluations = []

		for ind in pop:
			newInd = rankedIndividual.rankedIndividual(ind)
			newInd.setFitness(evalFitness(ind))
			fitnessEvaluations.append(newInd)

		sortedRanked = sorted(fitnessEvaluations , key = rankedIndividual.getFitness())

		for i in range(len(sortedRanked)):
			fitnessEvaluations[i].setProbability(float(i+1)/(totalRank))

		while(len(selected) < len(pop)):
			cummulativeSum = 0.0
			rand =  random.random()
			for j in len(fitnessEvaluations):
				cummulativeSum += fitnessEvaluations[i].getProbability()
				if(rand < cummulativeSum):
					selected.append(fitnessEvaluations[i].getIndividual)
					break

		return selected

	def ts(self , pop):
		selected = []

		while(len(selected) < len(pop)):
			ind1 = pop[random.randrange(len(pop))]
			ind2 = pop[random.randrange(len(pop))]

			f1 = evalFitness(ind1)
			f2 = evalFitness(ind1)

			if(f1 >f2):
				selected.append(ind1)
			else:
				selected.append(ind2)
		return selected
