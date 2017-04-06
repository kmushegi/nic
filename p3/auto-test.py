import os

def compile():
	cmd = "javac -Werror -Xlint:all ACO.java Node.java"
	os.system(cmd);

def runACS(numAnts, numIts, alpha, beta, rho, pfp, stopCond, secAllowed, 
			errAllowed, eps, qZero):
	cmd = "java ACO acs " + str(numAnts) + " " + str(numIts) + " " + str(alpha) \
			+ " " + str(beta) + " " + str(rho) + " " + fpf + " " + str(stopCond) \
			+ " " + str(secAllowed) + " " +str(errAllowed) + " " + str(eps) \
			+ " " + str(qZero) + " >> stats/acs-" + pfp + "-" + str(numAnts)  \
			+ "-" + str(numIts) + "-" + str(alpha) + "-" + str(beta) + "-" \
			+ "-" + str(rho) + "-" + str(eps) + "-" + str(qZero) + "-" \
			+ str(stopCond) + "-" + str(secAllowed) + "-" + str(errAllowed) \
			+ ".txt"

compile()


