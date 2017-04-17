JC = javac -Werror -Xlint:-serial
JV = java

ACORUNNER = ACORunner.java
ACO = ACO.java
UTILS = Utility.java
LOGGER = Logger.java
NODE = Node.java

#acs: 
#alg ants its alpha beta rho problemFP stopCond secAllowed errAllowed eps qZero

#eas:
#alg ants iters alpha beta rho problemFP stopCond secsAllowed errAllowed e

default: 
	$(JC) $(ACORUNNER) $(ACO) $(UTILS) $(LOGGER) $(NODE)

vis:
	$(JC) Visualizer.java

vis-run:
	@$(JV) Visualizer ../ALL_tsp/berlin52.tsp

acs:
	@$(JV) ACORunner acs 30 200 1 4 0.1 ../ALL_tsp/a280.tsp 0 30 0.2 0.1 0.9

eas:
	@$(JV) ACORunner eas 30 200 1 4 0.1 ../ALL_tsp/pr2392.tsp 0 15 0.1 30

clean:
	rm *.class