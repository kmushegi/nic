import os
import sys

problemdir = "sample-problems"

def compile():
	cmd = "javac -Werror -Xlint:all ACO.java Node.java"
	os.system(cmd)

def clear_stats_folder():
	cmd = "cd stats && rm *.txt"
	os.system(cmd)

def run(whichAlgorithm, numAnts, numIts, alpha, beta, rho, pfp, stopCond, secAllowed, 
			errAllowed, eps, qZero):
	cmd = "java ACO " + whichAlgorithm + " " + str(numAnts) + " " + str(numIts) \
			+ " " + str(alpha) + " " + str(beta) + " " + str(rho) + " " + pfp \
			+ " " + str(stopCond) + " " + str(secAllowed) + " " +str(errAllowed) \
			+ " " + str(eps) + " " + str(qZero) + " >> stats/acs-" \
			+ str(numAnts) + "-" + str(numIts) + "-" + str(alpha) + "-" + str(beta) \
			+ "-" + str(rho) + "-" + str(eps) + "-" + str(qZero) + "-" \
			+ str(stopCond) + "-" + str(secAllowed) + "-" + str(errAllowed) \
			+ ".txt"
	os.system(cmd)

compile()
clear_stats_folder()

numAnts = [5,10,20,30,40,50,60,80,100]
numIts = [100,500,1000,1500,2000,2500,3000]
alphas = [0,1,2,3,4,5]
betas = [0,1,2,3,4,5,6,7,8]
rhos = [0.0,0.1,0.2,0.3,0.4,0.5]
stopConditions = [0,1,2,3]
secondsAllowed = [10,20,30,40,50,60,70,80]
#erroAllowed = [...]
epsilons = [0.0,0.1,0.2,0.3,0.4,0.5]
qZeroes = [0.7,0.8,0.9,1.0,1.1,1.2]

whichAlgorithm = sys.argv[1]

for f in os.listdir(problemdir):
	if f.endswith(".tsp"):
		for n in numAnts:
			for i in numIts:
				for a in alphas:
					for b in betas:
						for r in rhos:
							for e in epsilons:
								for q in qZeroes:
									fp = problemdir + "/" + f
									run(whichAlgorithm,n,i,a,b,r,fp,0,10,0.1,e,q)









