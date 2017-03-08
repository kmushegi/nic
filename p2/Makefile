JC = javac -Werror -Xlint:all
JV = java

PSO = PSOTopologies.java
PARTICLE = Particle.java

default: 
	$(JC) $(PSO) $(PARTICLE)

run-default:
	@$(JV) PSOTopologies ri 49 rok

run-custom:
	@$(JV) PSOTopologies gl 49 5000 rok 25

clean:
	rm *.class
