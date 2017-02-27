JC = javac
JV = java

PSO = PSOTopologies.java

default: 
	$(JC) $(PSO)

run-default:
	@$(JV) PSOTopologies gl 50 rok

run-custom:
	@$(JV) PSOTopologies gl 50 5000 rok 25

clean:
	rm *.class
