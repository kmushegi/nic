"""
Neural Networks for Image Classification - Final Project
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Neural Networks for Image Classification, Final Project. This file
contains the implementation of a genetic algorithm with the ability implement rank selection 
and uniform crossovr. It also peforms the mutation operator. The genetic algorithm will find
the optimal set of parameters that maximizes accuracy of the neural network or the convolutonal 
neural network.
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

				#Unvaried parameters
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

			#For every parent pair in order of fitness after being ranked
			for i in range(0,len(parents),2):

				#Perform one point crossover
				if (self.crossover =="1c"):
					newChildren = self.onepoint(parents[i], parents[i+1])

				#Perform uniform crossover 
				elif (self.crossover =="uc"):
					newChildren = self.uniform(parents[i], parents[i+1])

				#Add the new children to the children container
				children.append(self.mutate(newChildren[0]))
				children.append(self.mutate(newChildren[1]))

			#Evaluate fitness of these children networks
			childFitnesses = []

			
			for child in children:

				#User specified neural network
				if (self.networkType == "nn"):
					childFitnesses.append(self.evaluateFitnessNN(child))

				#User specified convolutional neural network
				else:
					childFitnesses.append(self.evaluateFitnessCNN(child))

			#Find child netowork with the best fitness
			bestChild = children[childFitnesses.index(max(childFitnesses))]

			#Get actual fitness of this child network
			newBestFitness = max(childFitnesses)

			#Add this best fitness to "fitness over time"
			fitnessOT.append(newBestFitness)

			#Check if this is the new fittest neural net
			if (currentBestFitness < newBestFitness):

				#Keep track of this best neural network
				bestNN = bestChild

				#Keep track of this new best fitness
				currentBestFitness = newBestFitness

			#This group of children is now the new population
			self.population = children

			#Children fitness are now children for rank selection in next iteration
			nnFitnesses = childFitnesses

		#Print the best set of parameters to the stats file
		f = open(out_dest,'w')
		for p in bestNN.parameters:
			print >> f, p
		print >> f, fitnessOT
		f.close()

	#Rank selection
	def rs(self, nnFitnesses):
		#Willbe the selected networks
		selected = []

		#Will be used to calculate probabilites
		totalRank = ((self.populationSize)*(self.populationSize+1))/2

		#Will be a list of ranked individuals
		fitnessEvaluations = []

		#For every neural net 
		for i in range(len(self.population)):
			newInd = ri.rankedIndividual(self.population[i])
			newInd.setFitness(nnFitnesses[i])
			fitnessEvaluations.append(newInd)

		#Sort the indivuals by their fitness
		sorted(fitnessEvaluations, key=lambda rankedIndividual: rankedIndividual.fitness, reverse=True)

		#For every individual set its probability according to rank
		for i in range(len(fitnessEvaluations)):
			fitnessEvaluations[i].setProbability(float(i+1)/(totalRank))

		#Give every dingle individual a rank
		while(len(selected) < self.populationSize):

			#Sum of proabilities of individuals iterated through
			cummulativeSum = 0.0

			#Generate arandom number
			rand = random.random()

			#Get fitness of each individual
			for j in range(len(fitnessEvaluations)):

				#Add fitness of individual to cummulative sum
				cummulativeSum += fitnessEvaluations[j].getProbability()

				#Check if rand is less than cummulative sum so far
				if(rand < cummulativeSum):

					#If it is add this individual to selected pool
					selected.append(fitnessEvaluations[j].getIndividual())
					break

		return selected

	#Train and get accuracy which is fitness
	def evaluateFitnessNN(self, nn):
		return(nn.train(self.train_train_data[0], self.train_train_data[1])/1797.0)

	#Train and get accuracy which is fitness
	def evaluateFitnessCNN(self, cnn):
		return(cnn.train())

	def onepoint(self, p1, p2):
		return 1

	def uniform(self, p1, p2):

		#Perform crossover witth probability crossoverProb
		if (random.random() > self.crossoverProb):

			#If not the two children are equal to the two parents
			child1 = p1.parameters
			child2 = p2.parameters
			return (child1, child2)

		#Containers for the parameters that will initialize the networks
		child1 = []
		child2 = []

		#Loop through parameters in parents
		for key in range(len(p1.parameters)):

			#With probability 0.5 add to child one parameter from parent 1
			if (random.random() > 0.5):
				child1.append(p1.parameters[key])
			#Or else add the parameter from parent 2
			else:
				child1.append(p2.parameters[key])

		#Loop through the parameters in parents
		for key in range(len(p2.parameters)):

			#With probability 0.5 add to child two parameter from parent 1
			if (random.random() > 0.5):
				child2.append(p1.parameters[key])
			#Or else add the parameter from parent 2
			else:
				child2.append(p2.parameters[key])

		return (child1, child2)

	def mutate(self, child):

		#User specified neural network
		if (self.networkType == "nn"):

			#Child still set of parameters at this point so loop through parameters
			for i in range(len(child)):
				if (random.random() > self.mutationProb):

					#Unvaried parametera
					if (i == 0):
						n_in_neurons = 3072
						n_out_neurons = 10
						hidden_layer_info = random.choice(self.nnParams["hiddenInfo"])

						layer_info = [n_in_neurons]
						layer_info.extend(hidden_layer_info)
						layer_info.append(n_out_neurons)

						child[i] = layer_info

					#Choose a random choice of epoch parameters
					elif i == 1:
						child[i] = random.choice(self.nnParams["epochs"])
					#Choose a random chose of LR
					elif i == 2:
						child[i] = random.choice(self.nnParams["learningRate"])
					#Choose a random choice of starting weights
					else:
						child[i] = random.choice(self.nnParams["startWeights"])

			#Initialize the neural network
			net = nn.Network(child[0], child[1], child[2], child[3])

		#User specified a convolutional neural network
		else:
			#Childis still a set of parameters so loop through parameters
			for i in range(len(child)):
				if (random.random() > self.mutationProb):
					#Randomly choose each parameter for convolutional neural network
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

			#Initialize the convolutional neural network with these mutated parameters.
			net = cnn.CNetwork(child[0],child[1],child[2],child[3],child[4],child[5],
				child[6], self.x_y_train_data[0][0],self.x_y_train_data[0][1],
				self.x_y_train_data[1][0],self.x_y_train_data[1][1])
			net.build_network()

		return net
			


