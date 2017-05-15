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
		networkType = "nn"
		numIterations = 10
		populationSize = 20
		selection = "rs"
		crossover = "uc"
		crossoverProb = 0.8
		mutationProb = 0.2

	if (networkType == "nn"):
		nnParams = {
			"epochs": [10],
			"learningRate": [1, 0.5, 0.1],
			"hiddenInfo": [[], [100,50], []],
			"startWeights": [0,1,0.15, 0.2, 0.25],
		}

	else:
		nnParams = {
			"epochs": [1],
			"layers": [2,3,4],
			"dropout": [True, False],
			"batch_size": [32],
			"optimizer": ['sgd'],
			"data_augmentation": [True],
			"convActivation": ['relu'],
			"denseActivation": ['softmax']
		}

	geneticAlg = ga.GA(networkType, numIterations, populationSize, selection, crossover, crossoverProb, mutationProb, nnParams)
	geneticAlg.runGA()

main()




