import sys
import os

problemdir = "sample-problems"
statsdir = "stats/"

#compile the project
def compile():
	cmd = "javac EvolAlg.java"
	os.system(cmd)

#run Genetic Algorithm with given parameters and write statistics to an appropriately named file
def runGA(filename, numIndividuals, selection, crossover, crossprob, mutprob, numgenerations):
	cmd = "java EvolAlg " + filename + " " + str(numIndividuals) + " " + selection + " " + \
			str(crossover) + " " + str(crossprob) + " " + str(mutprob) + " " + str(numgenerations) + \
			" g >> " + "stats/ga-" + str(numIndividuals) + "-" + selection + "-" + crossover + "-" + \
			str(crossprob) + "-" + str(mutprob) + ".txt"
	os.system(cmd)
	os.system("echo \"\n\" >> " + "stats/ga-" + str(numIndividuals) + "-" + selection + "-" + crossover + \
			"-" + str(crossprob) + "-" + str(mutprob) + ".txt")

#run PBIL Algorithm with given parameters and write statistics to an appropriately named file
def runPBIL(filename, indsPerIter, posLR, negLR, mutprob, mutamt, numiterations):
	cmd = "java EvolAlg " + filename + " " + str(indsPerIter) + " " + str(posLR) + " " + \
			str(negLR) + " " + str(mutprob) + " " + str(mutamt) + " " + str(numiterations) + \
			" p >> " + "stats/pbil-" + str(indsPerIter) + "-" + str(posLR) + "-" + str(negLR) + \
			"-" + str(mutprob) + "-" + str(mutamt) + "-" + str(numiterations) + ".txt"
	os.system(cmd)
	os.system("echo \"\n\" >> " + "stats/pbil-" + str(indsPerIter) + "-" + str(posLR) + "-" + \
			str(negLR) + "-" + str(mutprob) + "-" + str(mutamt) + "-" + str(numiterations) + ".txt")

compile();

#GA Parameters to experiment with
selections = ["ts","rs","bs"]
crossovers = ["1c","uc"]
crossoverprobabilities = [0.6,0.7,0.8]
ga_mutationprobabilities = [0.01,0.05,0.1]

#GA Extended Params
cpe = [0.4,0.5,0.6,0.7,0.8]
mpe = [0.01,0.02,0.03,0.04,0.05,0.1]

#PBIL Parameters to experiment with
numberofindividuals = [100,175,250]
positivelearningrates = [0.1,0.25,0.4]
negativelearningrates = [0.075,0.25,0.4]
pbil_mutationprobabilities = [0.02,0.1,0.35]
mutationamounts = [0.05,0.2,0.4]

whichAlgorithm = sys.argv[1];

if(whichAlgorithm != "pbil" and whichAlgorithm != "ga" and \
	whichAlgorithm != "ga-extended" and whichAlgorithm != "pbil-extended" \
	and whichAlgorithm != "ga-extended-2"):
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
elif(whichAlgorithm == "ga-extended"): #ga-100-ts-uc-0.8-0.01.txt
	for c in cpe:
		for m in mpe:
			for filename in os.listdir(problemdir):
				if filename.endswith(".cnf"):
					fp = problemdir + "/" + filename
					runGA(fp,100,"ts","uc",c,m,1000)
				else:
					continue
elif(whichAlgorithm == "pbil-extended"): #pbil-250-0.1-0.075-0.1-0.05-1000.txt
	for filename in os.listdir(problemdir):
		if filename.endswith(".cnf"):
			fp = problemdir + "/" + filename
			runPBIL(fp,250,0.1,0.075,0.1,0.05,1000)
		else:
			continue
elif(whichAlgorithm == "ga-extended-2"): #ga-100-ts-uc-0.8-0.01.txt
	for filename in os.listdir(problemdir):
		if filename.endswith(".cnf"):
			fp = problemdir + "/" + filename
			runGA(fp,100,"ts","uc",0.8,0.01,1000)










