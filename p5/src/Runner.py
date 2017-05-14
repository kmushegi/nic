import ga as ga

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
		mutationProb = 0.01

	if (networkType == "nn"):
		nnParams = [
			"epochs": [],
			"layers": [],
			"startWeights": [],
			# "activation": [],
			"learningRate": []
		]

	else:
		nnParams = [
			"epochs": [],
			"layers": [],
			"startWeights": [],
			"activation": [],
			"learningRate": [],
			"dropout": [],
			"filters": [],
			"spatial": [],
			"stride": [],
			"padding":[]
		]


	geneticAlg = ga.GA(networkType, numIterations, populationSize, selection, crossover, crossoverProb, mutationProb, nnParams)
	geneticAlg.runGA()





