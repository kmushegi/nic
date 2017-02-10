JC = javac
JV = java

EA = EvolAlg.java

default: 
	$(JC) $(EA)

run:
	@$(JV) EvolAlg maxsat-problems/maxsat-random/highgirth/4SAT/HG-4SAT-V100-C900-4.cnf 100 0.1 0.075 0.02 0.05 1000 p

clean:
	rm *.class