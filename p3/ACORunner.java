/*
ACO for TSP - Project 3
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

The code in this file contains the implementation of Ant Colony Optimization 
for the Traveling Salesman Problem, as part of Project 3.
*/

import java.io.*;
import java.util.*;

public class ACORunner {

	public static void main(String[] args) {
		ACO acoInstance = new ACO();
		acoInstance.initializeACO(args);
		acoInstance.runACO();
	}
}