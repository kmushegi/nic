JC = javac -Werror -Xlint:-serial
JV = java

ACORUNNER = ACORunner.java
ACO = ACO.java
UTILS = Utility.java
LOGGER = Logger.java
NODE = Node.java
VIZ = Visualizer.java

#acs: 
#alg ants its alpha beta rho problemFP stopCond secAllowed errAllowed eps qZero vis(1/0)

#eas:
#alg ants iters alpha beta rho problemFP stopCond secsAllowed errAllowed e vis(1/0)

default: 
	$(JC) $(ACORUNNER) $(ACO) $(UTILS) $(LOGGER) $(VIZ) $(NODE)

acs:
	@$(JV) ACORunner acs 30 200 1 4 0.1 ../ALL_tsp/bier127.tsp 0 10 0.1 0.1 0.9 1

eas:
	@$(JV) ACORunner eas 30 200 1 4 0.1 ../ALL_tsp/pcb1173.tsp 0 15 0.1 30 1

clean:
	rm *.class