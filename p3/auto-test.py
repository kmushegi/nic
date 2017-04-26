"""
ACO for TSP - Project 3
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Ant Colony Optimization for the Traveling Salesman Problem,
Project 3. This file contains the implementation of our automated testing framework,
and is configured to run the ACS ant system.

Requires: joblib for parallelized execution.
"""

import os
import sys
from joblib import Parallel, delayed
import multiprocessing

problemdir = "sample-problems-2"

def compile():
	cmd = "cd src && make"
	os.system(cmd)

def clear_stats_folder():
	print "Cleaning stats folder"
	cmd = "cd stats && rm *.txt"
	os.system(cmd)

def runACS(whichAlgorithm, numAnts, numIts, alpha, beta, rho, pfp, stopCond, secAllowed, 
			errAllowed, eps, qZero):
	
	pfp = "../" + problemdir + "/" + pfp
	print pfp
	if not pfp.endswith(".tsp"):
		return
	else:
		cmd = "cd src && java ACORunner " + whichAlgorithm + " " + str(numAnts) + " " + str(numIts) \
			+ " " + str(alpha) + " " + str(beta) + " " + str(rho) + " " + pfp \
			+ " " + str(stopCond) + " " + str(secAllowed) + " " +str(errAllowed) \
			+ " " + str(eps) + " " + str(qZero) + " >> ../stats/" + whichAlgorithm + "-" \
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
parallel = int(sys.argv[2])

problems = os.listdir(problemdir)

if parallel:
	nj = multiprocessing.cpu_count()
	print "Spawning " + str(nj) + " jobs"
	Parallel(n_jobs=nj)(delayed(run)(whichAlgorithm=whichAlgorithm,numAnts=ni,numIts=ii,alpha=ai,beta=bi,rho=ri,pfp=fpi,stopCond=0,secAllowed=10,errAllowed=0.1,eps=ei,qZero=qi) for ni in numAnts for ii in numIts for ai in alphas for bi in betas for ri in rhos for ei in epsilons for qi in qZeroes for fpi in problems)
	print "Done"
else:
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
										runACS(whichAlgorithm,n,i,a,b,r,fp,0,10,0.1,e,q)
									










