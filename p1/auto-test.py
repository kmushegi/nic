import sys
import os

problemdir = "sample-problems"
statisdir = "stats/"

def compile():
	cmd = "javac EvolAlg.java"
	os.system(cmd)

def runGA(filename, numIndividuals, selection, crossover, crossprob, mutprob, numgenerations):
	cmd = "java EvolAlg " + filename + " " + str(numIndividuals) + " " + selection + " " + \
			str(crossover) + " " + str(crossprob) + " " + str(mutprob) + " " + str(numgenerations) + \
			" g >> " + "stats/ga-" + str(numIndividuals) + "-" + selection + "-" + crossover + "-" + \
			str(crossprob) + "-" + str(mutprob) + ".txt"
	os.system(cmd)
	os.system("echo \"\\n\" >> " + "stats/ga-" + str(numIndividuals) + "-" + selection + "-" + crossover + \
			"-" + str(crossprob) + "-" + str(mutprob) + ".txt")

def runPBIL(filename, indsPerIter, posLR, negLR, mutprob, mutamt, numiterations):
	cmd = "java EvolAlg " + filename + " " + str(indsPerIter) + " " + str(posLR) + " " + \
			str(negLR) + " " + str(mutprob) + " " + str(mutamt) + " " + str(numiterations) + \
			" p >> " + "stats/pbil-" + str(indsPerIter) + "-" + str(posLR) + "-" + str(negLR) + \
			"-" + str(mutprob) + "-" + str(mutamt) + "-" + str(numiterations) + ".txt"
	os.system(cmd)
	os.system("echo \"\n\" >> " + "stats/pbil-" + str(indsPerIter) + "-" + str(posLR) + "-" + \
			str(negLR) + "-" + str(mutprob) + "-" + str(mutamt) + "-" + str(numiterations) + ".txt")

compile();

#GA Parameters
selections = ["ts","rs","bs"]
crossovers = ["1c","uc"]
crossoverprobabilities = [0.5,0.6,0.7,0.8,0.9]
ga_mutationprobabilities = [0.0025,0.005,0.01,0.05,0.1,0.15]

#PBIL Parameters
numberofindividuals = [50,100,150,200,250]
positivelearningrates = [0.025,0.05,0.1,0.15,0.2]
negativelearningrates = [0.025,0.05,0.075,0.1,0.125]
pbil_mutationprobabilities = [0.01,0.02,0.03,0.04,0.05]
mutationamounts = [0.03,0.04,0.05,0.06,0.07]

whichAlgorithm = sys.argv[1];

if(whichAlgorithm == "pbil"):
	print "Running PBIL"
elif(whichAlgorithm == "ga"):
	print "Running GA"
else:
	print "Argument Specified Incorrectly. Supply 'pbil' or 'ga'"
	sys.exit()


if(whichAlgorithm == "pbil"):
	for n in numberofindividuals:
		for pl in positivelearningrates:
			for nl in negativelearningrates:
				for mutp in pbil_mutationprobabilities:
					for mtamt in mutationamounts:
						for filename in os.listdir(problemdir):
							if filename.endswith(".cnf"):
								fp = problemdir + "/" + filename
								runPBIL(fp,n,pl,nl,mutp,mtamt,1000)
							else:
								continue
elif(whichAlgorithm == "ga"):
	for s in selections:
		for crs in crossovers:
			for cp in crossoverprobabilities:
				for mutp in ga_mutationprobabilities:
					for filename in os.listdir(problemdir):
						if filename.endswith(".cnf"):
							fp = problemdir + "/" + filename
							runGA(fp,100,s,crs,cp,mutp,1000)
						else:
							continue








