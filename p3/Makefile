JC = javac -Werror -Xlint:all
JV = java

ACO = ACO.java
NODE = Node.java

#acs: 
#alg ants its alpha beta rho problemFP stopCond secAllowed errAllowed eps qZero

#eas:
#alg ants iters alpha beta rho problemFP stopCond secsAllowed errAllowed e

default: 
	$(JC) $(ACO) $(NODE)

acs:
	@$(JV) ACO acs 30 200 1 4 0.1 ALL_tsp/a280.tsp 0 10 0.1 0.1 0.9

eas:
	@$(JV) ACO eas 30 200 1 3 0.1 ALL_tsp/berlin52.tsp 0 15 0.1 30

clean:
	rm *.class