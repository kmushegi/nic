JC = javac
JV = java

EA = EvolAlg.java

default: 
	$(JC) $(EA)

p:
	@$(JV) EvolAlg maxsat-problems/maxsat-random/max3sat/60v/s3v60c800-1.cnf 100 0.1 0.075 0.02 0.05 1000 p

g:
	@$(JV) EvolAlg sample-problems/HG-3SAT-V250-C1000-1.cnf 100 bs 1c 0.6 0.01 1000 g

clean:
	rm *.class
