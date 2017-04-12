/*
ACO for TSP - Project 3
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Ant Colony Optimization for the Traveling Salesman Problem,
Project 3. This file is the point of entry for our ACO implementation. It is
responsible for creating an instance of an ACO system and running it using
the parameters supplied on the command line.
*/

public class ACORunner {

	public static void main(String[] args) {
		ACO acoInstance = new ACO();
		acoInstance.initializeACO(args);
		acoInstance.runACO();
	}
}