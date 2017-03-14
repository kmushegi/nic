JC = javac -Werror -Xlint:all
JV = java

PSO = PSOTopologies.java
PARTICLE = Particle.java

default: 
	$(JC) $(PSO) $(PARTICLE)

run-default:

	@$(JV) PSOTopologies vn 50 ack

run-custom:
	@$(JV) PSOTopologies ra 49 100000 rok 10

clean:
	rm *.class
