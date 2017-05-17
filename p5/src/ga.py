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

out_dest = '../stats/stats-1.txt'

class GA(object):

	#Constructor
	def __init__(self, networkType, numIterations, populationSize, selection, crossover, crossoverProb, mutationProb, nnParams):
		self.networkType = networkType
		self.numberOfIterations = numIterations
		self.populationSize = populationSize
		self.selection = selection
		self.crossover = crossover
		self.crossoverProb = crossoverProb
		self.mutationProb = mutationProb
		self.nnParams = nnParams

		#User specified neural network
		if (self.networkType == "nn"):
			sys.stdout.flush()

			#Load data from cifar10 dataset and 10 output neurons
			self.train_train_data = dl.get_data("nn",'cifar10',10)

		#User specified convolutional neural network
		else:
			sys.stdout.flush()

			#Load data from cifar10
			self.x_y_train_data = dl.get_data("cnn",'cifar10')

		self.population = self.initializePopulation()

	def initializePopulation(self):

		#Container for the networks
		pop = []
		for _ in range(0, self.populationSize):

			#User specified neural network
			if (self.networkType == "nn"):

				#Number of input nodes from rgb values
				n_in_neurons = 3072

				#Number of output nodes for 10 categories
				n_out_neurons = 10

				#Randomly hoose hidden layer parameters
				hidden_layer_info = random.choice(self.nnParams["hiddenInfo"])


				layer_info = [n_in_neurons]
				layer_info.extend(hidden_layer_info)
				layer_info.append(n_out_neurons)

				#Crete the neural net with specified node info and random parameters; specified in network.py
				net = nn.Network(layer_info, random.choice(self.nnParams["epochs"]),
				random.choice(self.nnParams["learningRate"]), random.choice(self.nnParams["startWeights"])) #Parameters

			#User specified convolutional neural network
			else:

				#Initialize convolution neural net with random parameters; specified in cnetwork.py
				net = cnn.CNetwork(random.choice(self.nnParams["epochs"]), random.choice(self.nnParams["dropout"]),
				random.choice(self.nnParams["batch_size"]), random.choice(self.nnParams["optimizer"]), random.choice(self.nnParams["data_augmentation"]),
				random.choice(self.nnParams["convActivation"]), random.choice(self.nnParams["denseActivation"]),
				self.x_y_train_data[0][0],self.x_y_train_data[0][1],self.x_y_train_data[1][0],self.x_y_train_data[1][1]) #Parameters
				net.build_network()

			#Append initialized network to population container
			pop.append(net)

		#Return randomly initialized neural network population
		return pop

	def runGA(self):

		#Container for the fitness(accuracy) of each network
		nnFitnesses = []
		for i in range(len(self.population)):

			#User specified neural network
			if (self.networkType == "nn"):
				nnFitnesses.append(self.evaluateFitnessNN(self.population[i]))

			#User specified convolutional neural network
			else:
				nnFitnesses.append(self.evaluateFitnessCNN(self.population[i]))

		#Find network with best fitness
		bestNN = self.population[nnFitnesses.index(max(nnFitnesses))]

		#Retrieve the actual fitness of this network
		currentBestFitness = max(nnFitnesses)

		#Keeps track of the best fitness over number of iterations
		fitnessOT = []
		fitnessOT.append(currentBestFitness)

		#Container for networks selected for crossover
		parents = []

		#For number of iterations
		for i in range(self.numberOfIterations):

			#Prints iteration number to stats file
			f = open(out_dest,'w')
			print >> f, i
			f.close()

			#If selected rank selection run rank selection
			if self.selection == "rs":
				parents = self.rs(nnFitnesses)
			else:
				sys.exit(1)

			#Container for the children that resulting from crossover
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
		f.close()

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
						child[i] = random.choice(self.nnParams["dropout"])
					elif i == 2:
						child[i] = random.choice(self.nnParams["batch_size"])
					elif i == 3:
						child[i] = random.choice(self.nnParams["optimizer"])
					elif i == 4:
						child[i] = random.choice(self.nnParams["data_augmentation"])
					elif i == 5:
						child[i] = random.choice(self.nnParams["convActivation"])
					elif i == 6:
						child[i] = random.choice(self.nnParams["denseActivation"])

			net = cnn.CNetwork(child[0],child[1],child[2],child[3],child[4],child[5],
				child[6], self.x_y_train_data[0][0],self.x_y_train_data[0][1],
				self.x_y_train_data[1][0],self.x_y_train_data[1][1])
			net.build_network()

		return net
			


