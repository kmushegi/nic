import os
import sys

# from joblib import Parallel, delayed
# import multiprocessing

problemdir = "sample-problems"

def compile():
	cmd = "cd src && make"
	os.system(cmd)

def clear_stats_folder():
	print "Cleaning stats folder"
	cmd = "cd stats && rm *.txt"
	os.system(cmd)

def run(whichAlgorithm, numAnts, numIts, alpha, beta, rho, pfp, stopCond, secAllowed, 
			errAllowed, eps, qZero):
	
	if not pfp.endswith(".tsp"):
		return
	else:
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

numAnts = [10]

whichAlgorithm = sys.argv[1]
# parallel = int(sys.argv[2])

problems = os.listdir(problemdir)

# if parallel:
# 	nj = multiprocessing.cpu_count()
# 	print "Spawning " + str(nj) + " jobs"
# 	Parallel(n_jobs=nj)(delayed(run)(whichAlgorithm=whichAlgorithm,numAnts=ni,numIts=ii,alpha=ai,beta=bi,rho=ri,pfp=fpi,stopCond=0,secAllowed=10,errAllowed=0.1,eps=ei,qZero=qi) for ni in numAnts for ii in numIts for ai in alphas for bi in betas for ri in rhos for ei in epsilons for qi in qZeroes for fpi in problems)
# 	print "Done"
# else:
for i in range(0,3):
	for n in numAnts:
		for f in os.listdir(problemdir):
			if f.endswith(".tsp"):
				fp = "../" + problemdir + "/" + f
				run(whichAlgorithm,n,200,1,4,0.1,fp,0,10,0.1,0.1,0.9)
									










