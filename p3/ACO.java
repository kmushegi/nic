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

public class ACO {

	//*******************************************
	//*Parameters in order of acceptance from CL*
	//*******************************************

	//General parameters
	private static int whichAlgorithm; // 1 - acs; 0 - eas
	private static int numberofAnts;
	private static int numberofIterations;
	private static double alpha;
	private static double beta;
	private static double rho;
	private static String problemFilePath; //location of .tsp file

	//ACS Parameters
	private static double eps;
	private static double tauZero; //calculated not provided
	private static double qZero;

	/* τ0 = 1/(n·Lnn) (ACS) Lnn is the length of a nearest neighbor tour,
	 a tour constructed by starting at an arbitrary city and selecting 
	 the closest unvisited city to visit next (note that nearest neighbor 
	 tours are not unique) 
	*/

	//EAS Parameters
	private static double e;

	//ACO Constants
	private static final String acs = "acs";
	private static final String eas = "eas";
	private static final String incorrectParams = "Parameters were not supplied correctly";

	public static void main(String[] args) {
		readParams(args);
		printParams(whichAlgorithm);
	}

	public static void readParams(String[] args) {
		if(args.length != 9 && args.length != 8) {
			printErrorAndExit();
		} 

		if(args[0].equals(acs)) {
			whichAlgorithm = 1;
			eps = Double.parseDouble(args[6]);
			qZero = Double.parseDouble(args[7]);
			problemFilePath = args[8];
		} else if(args[0].equals(eas)) {
			whichAlgorithm = 0;
			e = Double.parseDouble(args[6]);
			problemFilePath = args[7];
		} else {
			printErrorAndExit();
		}

		numberofAnts = Integer.parseInt(args[1]);
		numberofIterations = Integer.parseInt(args[2]);
		alpha = Double.parseDouble(args[3]);
		beta = Double.parseDouble(args[4]);
		rho = Double.parseDouble(args[5]);
	}

	public static void printParams(int alg) {
		System.out.println("Algorithm: " + ((alg == 1) ? acs : eas)
			+ "\n# of Ants: " + numberofAnts
			+ "\n# of Iterations: " + numberofIterations
			+ "\nAlpha: " + alpha
			+ "\nBeta: " + beta
			+ "\nRho: " + rho);

		if(alg == 1) { //print acs params
			System.out.println("Epsilon: " + eps
				+ "\ntauZero: " + tauZero
				+ "\nqZero: " + qZero);
		} else if(alg == 0) { //print eas params
			System.out.println("E: " + e);
		}
		System.out.println("Problem File: " + problemFilePath);
	}

	public static void printErrorAndExit() {
		System.out.println(incorrectParams);
		System.exit(1); //exit with error
	}
}