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

import java.util.*;

public class ACORunner {

	public static void main(String[] args) {
		ACO acoInstance = new ACO(); //create ACO instance
		acoInstance.initializeACO(args); //initialize the object with given args
		ArrayList<ArrayList<Node>> s = acoInstance.runACO(); //run the algorithm with given args

		if(Integer.parseInt(args[args.length - 1]) == 1) { //if visualization on, visualize
			Visualizer vis = new Visualizer();
			vis.display(s.get(0),s.get(1));
		}
	}
}