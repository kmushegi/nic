JC = javac
JV = java

PSO = PSOTopologies.java
PARTICLE = Particle.java

default: 
	$(JC) $(PSO) $(PARTICLE)

run-default:
	@$(JV) PSOTopologies gl 50 rok

run-custom:
	@$(JV) PSOTopologies gl 50 5000 rok 25

clean:
	rm *.class
