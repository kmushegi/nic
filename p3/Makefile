JC = javac -Werror -Xlint:all
JV = java

ACO = ACO.java

default: 
	$(JC) $(ACO)

acs:
	@$(JV) ACO acs 30 1000 1 3 0.1 0.1 0.9

eas:
	@$(JV) ACO eas 30 1000 1 3 0.1 30

clean:
	rm *.class