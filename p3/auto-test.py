import os
import sys

problemdir = "sample-problems"

def compile():
	cmd = "cd src && make"
	os.system(cmd)

def clear_stats_folder():
	cmd = "cd stats && rm *.txt"
	os.system(cmd)

def run(whichAlgorithm, numAnts, numIts, alpha, beta, rho, pfp, stopCond, secAllowed, 
			errAllowed, eps, qZero):
	cmd = "cd src && java ACORunner " + whichAlgorithm + " " + str(numAnts) + " " + str(numIts) \
			+ " " + str(alpha) + " " + str(beta) + " " + str(rho) + " " + pfp \
			+ " " + str(stopCond) + " " + str(secAllowed) + " " +str(errAllowed) \
			+ " " + str(eps) + " " + str(qZero) + " >> ../stats/acs-" \
			+ str(numAnts) + "-" + str(numIts) + "-" + str(alpha) + "-" + str(beta) \
			+ "-" + str(rho) + "-" + str(eps) + "-" + str(qZero) + "-" \
			+ str(stopCond) + "-" + str(secAllowed) + "-" + str(errAllowed) \
			+ ".txt"
	os.system(cmd)

compile()
clear_stats_folder()

numAnts = [20,30,40]
numIts = [100,200,300]
alphas = [1,2,3]
betas = [2,3,4]
rhos = [0.1,0.2,0.3]
epsilons = [0.1,0.2,0.3]
qZeroes = [0.7,0.8,0.9]

whichAlgorithm = sys.argv[1]

for n in numAnts:
	for i in numIts:
		for a in alphas:
			for b in betas:
				for r in rhos:
					for e in epsilons:
						for q in qZeroes:
							for f in os.listdir(problemdir):
								if f.endswith(".tsp"):
									fp = "../" + problemdir + "/" + f
									run(whichAlgorithm,n,i,a,b,r,fp,0,10,0.1,e,q)









