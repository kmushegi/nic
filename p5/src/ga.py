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
import rankedIndividual as ri
import data_loader as dl

out_dest = '/home/kmushegi/p5/stats/stats.txt'

class GA(object):

	def __init__(self, networkType, numIterations, populationSize, selection, crossover, crossoverProb, mutationProb, nnParams):
		self.networkType = networkType
		self.numberOfIterations = numIterations
		self.populationSize = populationSize
		self.selection = selection
		self.crossover = crossover
		self.crossoverProb = crossoverProb
		self.mutationProb = mutationProb
		self.nnParams = nnParams

		if (self.networkType == "nn"):
			sys.stdout.flush()
			self.train_train_data = dl.get_data("nn",'cifar10',10) #Length = 2

		else:
			sys.stdout.flush()
			self.x_y_train_data = dl.get_data("cnn",'cifar10') #Array of arrays

		self.population = self.initializePopulation()

	def initializePopulation(self):
		pop = []
		for _ in range(0, self.populationSize):
			if (self.networkType == "nn"):
				n_in_neurons = 3072
				n_out_neurons = 10
				hidden_layer_info = random.choice(self.nnParams["hiddenInfo"])

				layer_info = [n_in_neurons]
				layer_info.extend(hidden_layer_info)
				layer_info.append(n_out_neurons)

				net = nn.Network(layer_info, random.choice(self.nnParams["epochs"]),
				random.choice(self.nnParams["learningRate"]), random.choice(self.nnParams["startWeights"])) #Parameters

			else:
				net = cnn.CNetwork(random.choice(self.nnParams["epochs"]), random.choice(self.nnParams["layers"]), random.choice(self.nnParams["dropout"]),
				random.choice(self.nnParams["batch_size"]), random.choice(self.nnParams["optimizer"]), random.choice(self.nnParams["data_augmentation"]),
				random.choice(self.nnParams["convActivation"]), random.choice(self.nnParams["denseActivation"]),
				self.x_y_train_data[0][0],self.x_y_train_data[0][1],self.x_y_train_data[1][0],self.x_y_train_data[1][1]) #Parameters
				net.build_network()

			pop.append(net)

		return pop

	def runGA(self):

		nnFitnesses = []
		for i in range(len(self.population)):
			if (self.networkType == "nn"):
				nnFitnesses.append(self.evaluateFitnessNN(self.population[i]))
			else:
				nnFitnesses.append(self.evaluateFitnessCNN(self.population[i]))

		bestNN = self.population[nnFitnesses.index(max(nnFitnesses))]
		currentBestFitness = max(nnFitnesses)

		fitnessOT = []
		fitnessOT.append(currentBestFitness)

		parents = []

		for _ in range(self.numberOfIterations):
			if self.selection == "rs":
				parents = self.rs(nnFitnesses)
			else:
				sys.exit(1)

			children = []

			for i in range(0,len(parents),2):
				if (self.crossover =="1c"):
					newChildren = self.onepoint(parents[i], parents[i+1])
				elif (self.crossover =="uc"):
					newChildren = self.uniform(parents[i], parents[i+1])
				
				children.append(self.mutate(newChildren[0]))
				children.append(self.mutate(newChildren[1]))

			childFitnesses = []

			for child in children:
				if (self.networkType == "nn"):
					childFitnesses.append(self.evaluateFitnessNN(child))
				else:
					childFitnesses.append(self.evaluateFitnessCNN(child))

			bestChild = children[childFitnesses.index(max(childFitnesses))]
			newBestFitness = max(childFitnesses)

			fitnessOT.append(newBestFitness)

			if (currentBestFitness < newBestFitness):
				bestNN = bestChild
				currentBestFitness = newBestFitness

			self.population = children
			nnFitnesses = childFitnesses


		f = open(out_dest,'w')
		for p in bestNN.parameters:
			print >> f, p
		print >> f, fitnessOT

	def rs(self, nnFitnesses):
		#Willbe the selected networks
		selected = []

		#Willbe used to calculate probabilites
		totalRank = ((self.populationSize)*(self.populationSize+1))/2

		#Will be a list of ranked individuals
		fitnessEvaluations = []

		for i in range(len(self.population)):
			newInd = ri.rankedIndividual(self.population[i])
			newInd.setFitness(nnFitnesses[i])
			fitnessEvaluations.append(newInd)

		sorted(fitnessEvaluations, key=lambda rankedIndividual: rankedIndividual.fitness, reverse=True)

		for i in range(len(fitnessEvaluations)):
			fitnessEvaluations[i].setProbability(float(i+1)/(totalRank))

		while(len(selected) < self.populationSize):
			cummulativeSum = 0.0
			rand = random.random()
			for j in range(len(fitnessEvaluations)):
				cummulativeSum += fitnessEvaluations[j].getProbability()
				if(rand < cummulativeSum):
					selected.append(fitnessEvaluations[j].getIndividual())
					break

		return selected

	def evaluateFitnessNN(self, nn):
		return(nn.train(self.train_train_data[0], self.train_train_data[1])/1797.0)

	def evaluateFitnessCNN(self, cnn):
		return(cnn.train())

	def onepoint(self, p1, p2):
		return 1

	def uniform(self, p1, p2):
		if (random.random() > self.crossoverProb):
			child1 = p1.parameters
			child2 = p2.parameters
			return (child1, child2)

		child1 = []
		child2 = []

		for key in range(len(p1.parameters)):
			if (random.random() > 0.5):
				child1.append(p1.parameters[key])
			else:
				child1.append(p2.parameters[key])

		for key in range(len(p2.parameters)):
			if (random.random() > 0.5):
				child2.append(p1.parameters[key])
			else:
				child2.append(p2.parameters[key])

		return (child1, child2)

	def mutate(self, child):

		if (self.networkType == "nn"):
			for i in range(len(child)):
				if (random.random() > self.mutationProb):

					if (i == 0):
						n_in_neurons = 3072
						n_out_neurons = 10
						hidden_layer_info = random.choice(self.nnParams["hiddenInfo"])

						layer_info = [n_in_neurons]
						layer_info.extend(hidden_layer_info)
						layer_info.append(n_out_neurons)

						child[i] = layer_info
					elif i == 1:
						child[i] = random.choice(self.nnParams["epochs"])
					elif i == 2:
						child[i] = random.choice(self.nnParams["learningRate"])
					else:
						child[i] = random.choice(self.nnParams["startWeights"])

			net = nn.Network(child[0], child[1], child[2], child[3])

		else:
			for i in range(len(child)):
				if (random.random() > self.mutationProb):
					if (i == 0):
						child[i] = random.choice(self.nnParams["epochs"])
					elif i == 1:
						child[i] = random.choice(self.nnParams["layers"])
					elif i == 2:
						child[i] = random.choice(self.nnParams["dropout"])
					elif i == 3:
						child[i] = random.choice(self.nnParams["batch_size"])
					elif i == 4:
						child[i] = random.choice(self.nnParams["optimizer"])
					elif i == 5:
						child[i] = random.choice(self.nnParams["data_augmentation"])
					elif i == 6:
						child[i] = random.choice(self.nnParams["convActivation"])
					elif i == 7:
						child[i] = random.choice(self.nnParams["denseActivation"])

			net = cnn.CNetwork(child[0],child[1],child[2],child[3],child[4],child[5],
				child[6],child[7], self.x_y_train_data[0][0],self.x_y_train_data[0][1],
				self.x_y_train_data[1][0],self.x_y_train_data[1][1])
			net.build_network()

		return net
			


