JC = javac
JV = java

EA = EvolAlg.java

default: 
	$(JC) $(EA)

p:
	@$(JV) EvolAlg maxsat-problems/maxsat-random/max3sat/60v/s3v60c800-1.cnf 100 0.1 0.075 0.02 0.05 1000 p

g:
	@$(JV) EvolAlg maxsat-problems/maxsat-random/max3sat/60v/s3v60c800-1.cnf 100 bs uc 0.7 0.01 1000 g

clean:
	rm *.class