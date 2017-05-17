"""
Evolving Neural Networks Using Genetic Algorithms - Project 5
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Project 5. This file is an entry point into the program and let's the user
run the Genetic Algorithm to evolve optimal parameters for the Neural Networks (MLP, CNN).

If the Command Line parameters are correctly specified (lines 22-28) then the parameters are read
from the CL, otherwise Genetic Algorithm parameters default to ones defined on lines 32-38.
"""

import ga as ga
import data_loader as dl
import sys
import network as nn

def main():
	if len(sys.argv) == 6:
		networkType = sys.argv[1]
		numIterations = int(sys.argv[2])
		populationSize = int(sys.argv[3])
		selection = sys.argv[4]
		crossover = sys.argv[5]
		crossoverProb = float(sys.argv[6])
		mutationProb = float(sys.argv[7])

	else:
		networkType = "cnn"
		numIterations = 1
		populationSize = 2
		selection = "rs"
		crossover = "uc"
		crossoverProb = 0.8
		mutationProb = 0.2

	if (networkType == "nn"):
		nnParams = {
			"epochs": [10],
			"learningRate": [1, 0.5, 0.1],
			"hiddenInfo": [[], [300], [150,100]],
			"startWeights": [0,1,0.15, 0.2, 0.25],
		}

	else:
		nnParams = {
			"epochs": [1],
			"dropout": [True, False],
			"batch_size": [32],
			"optimizer": ['sgd','rmsprop','adam'],
			"data_augmentation": [False],
			"convActivation": ['relu','elu','sigmoid'],
			"denseActivation": ['softmax','tanh','sigmoid']
		}

	geneticAlg = ga.GA(networkType, numIterations, populationSize, selection, crossover, crossoverProb, mutationProb, nnParams)
	geneticAlg.runGA()

main()




