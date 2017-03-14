import sys
import os

def compile():
	cmd = "javac -Werror -Xlint:all PSOTopologies.java Particle.java"
	os.system(cmd)

def runPSO(whichTopology, swarmSize, iterations, whichFunction, dimensions):
	cmd = "java PSOTopologies " + whichTopology + " " + str(swarmSize) + " " + \
			str(iterations) + " " + whichFunction + " " + str(dimensions) + \
			" >> stats/pso-" + whichTopology + "-" + str(swarmSize) + \
			"-" + str(iterations) + "-" + whichFunction + "-" + str(dimensions) + ".txt"
	os.system(cmd)


compile()

swarm_size = [16,30,49]
topology = ["ri","vn","ra","gl"]
functions = ["rok","ack","ras"]

dimensions = 30
iterations = 10000
runs = 20

for s in swarm_size:
	for t in topology:
		for f in functions:
			for r in range(runs):
				runPSO(t,s,iterations,f,dimensions)
