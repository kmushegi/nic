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

	private static int stopCondition;
	// 0 - terminate after max iterations reached
	// 1 - terminate if found tour that is no more than a specified percentage 
	//	over the optimal (0.0 would mean you will not settle for anything less than the optimal)
	// 2 - both
	// 3 - terminate after secondsAllowed exceeds
	// 4 - all three

	private static int secondsAllowed ;
	private static int optimalTourCost;
	private static double errorAllowed;

	//ACS parameters
	private static double eps;
	private static double tauZero; //calculated not provided
	private static double qZero;

	/* τ0 = 1/(n·Lnn) (ACS) Lnn is the length of a nearest neighbor tour,
	 a tour constructed by starting at an arbitrary city and selecting 
	 the closest unvisited city to visit next (note that nearest neighbor 
	 tours are not unique) TODO;
	*/

	//EAS parameters
	private static double e;

	//*******************************************

	//Containers
	private static ArrayList<Node> nodes = new ArrayList<>();
	private static ArrayList<Node> solutionTour;
	private static int solutionCost;
	private static double[][] pheromoneMatrix;

	//line below and all new's up here will be moved to the constructor of ACO class
	private static Map<String, Integer> optTourLengths = new HashMap<String, Integer>();

	//ACO constants
	private static final String acs = "acs";
	private static final String eas = "eas";
	private static final String optTourLengthsFilePath = "optimalTourLengths.txt";
	private static final String incorrectParams = "Parameters were not supplied correctly";

	//Util
	private static Random generator = new Random();

	public static void main(String[] args) {
		readParams(args);
		readProblem(problemFilePath);
		readOptTourLengths(optTourLengthsFilePath);
		printParams(whichAlgorithm);

		if(whichAlgorithm == 1) {
			solutionTour = acs(numberofAnts, numberofIterations, alpha, beta, 
								rho, eps, tauZero, qZero);
		} else if(whichAlgorithm == 0) {
			solutionTour = eas(numberofAnts, numberofIterations, alpha, beta, 
								rho, e);
		} else {
			printErrorAndExit();
		}
		solutionCost = computeCost(solutionTour);
		outputSolution(solutionTour, solutionCost);
	}

	public static ArrayList<Node> acs(int ants, int its, double a, double b,
									double r, double eps, double t, double q) {
		ArrayList<Node> bestTour = initializeRandomSolution(nodes);
		int bestCost = computeCost(bestTour);
		double initialPheromone = 1.0 / ((double)nodes.size() * bestCost);
		pheromoneMatrix = initializePheromoneMatrix(nodes.size(), initialPheromone);

		int itCounter = 0;
		long st = System.nanoTime();

		while(evaluateStopCondition(stopCondition,itCounter, st, bestTour)) {
			itCounter += 1;

			//store solutions here?
			for(int i = 0; i < nodes.size(); i++) {
				ArrayList<Node> candidateTour;
				int candidateCost;

				// candidateTour = constructSolution(pheromoneMatrix, b, q);

			}
		}

		return bestTour;
	}

	public static ArrayList<Node> eas(int ants, int its, double a, double b,
									double r, double e) {
		ArrayList<Node> bestTour = new ArrayList<>();

		//do stuff

		return bestTour;
	}

	public static ArrayList<Node> constructSolution(double[][] pm, double b, double q) {
		ArrayList<Integer> tour =  new ArrayList<>();
		tour.add(generator.nextInt(nodes.size() + 1));

		while(tour.size() != nodes.size()) {

		}

		ArrayList<Node> tmp = new ArrayList<>();
		return tmp;
	}


	public static Boolean evaluateStopCondition (int stopCondition, int currIt, 
										long startTime, ArrayList<Node> tour) {

		switch(stopCondition) {
			case 0: 
				return (currIt < numberofIterations);
			case 1:	
				return (((computeCost(tour) / optimalTourCost) - 1) > errorAllowed);
			case 2:
				return ((currIt < numberofIterations) 
					&& (((computeCost(tour) / optimalTourCost) - 1) > errorAllowed));
			case 3:
				return (((System.nanoTime() - startTime) / 1000000000.0) < secondsAllowed);
			case 4:
				return ((currIt < numberofIterations) 
					&& (((computeCost(tour) / optimalTourCost) - 1) > errorAllowed) 
					&& ((System.nanoTime() - startTime) / 1000000000.0) < secondsAllowed);
			default: printErrorAndExit();
		}
		return false;
	}

	// 0 - terminate after max iterations reached
	// 1 - terminate if found tour that is no more than a specified percentage 
	//	over the optimal (0.0 would mean you will not settle for anything less than the optimal)
	// 2 - both
	// 3 - terminate after secondsAllowed exceeds
	// 4 - all three

	public static double euclideanDistance2D(Node n1, Node n2) {
		return Math.sqrt(Math.pow(n1.x-n2.x,2) + Math.pow(n1.y-n2.y,2));
	}

	public static int computeCost(ArrayList<Node> t) {
		int distance = 0;
		for(int i = 0; i < t.size(); i++) {
			Node n1 = t.get(i);
			Node n2 = ((i == t.size() - 1) ? t.get(0) : t.get(i+1));
			distance += (euclideanDistance2D(n1,n2));
		}
		return distance;
	}

	public static double[][] initializePheromoneMatrix(int n, double initPh) {
		double[][] temp = new double[n][n];

		for(int r = 0; r < n; r++) {
			for(int c = 0; c < n; c++) {
				temp[r][c] = initPh;
			}
		}
		return temp;
	}

	public static ArrayList<Node> initializeRandomSolution(ArrayList<Node> s) {
		ArrayList<Node> randSol = s;
		long seed = System.nanoTime();
		Collections.shuffle(randSol, new Random(seed));
		return randSol;
	}

	public static void processProblemLine(String[] tokens) {
		try {
			Integer.parseInt(tokens[0]); //make sure its a node line
			tokens = formatNodeInput(tokens);

			Node temp = new Node(Integer.parseInt(tokens[0]),
							Integer.parseInt(tokens[1]),
							Integer.parseInt(tokens[2]));
			nodes.add(temp);
		} catch (NumberFormatException e) {
			// System.err.format("NumberFormatException %s%n",e);
		}
	}

	public static void processOptTourLine(String[] tokens) {
		try {
			optTourLengths.put(tokens[0],Integer.parseInt(tokens[1]));
		} catch (NumberFormatException e) {
			// System.err.format("NumberFormatException %s%n",e);
		}
	}

	public static void readProblem(String fp) {
		try(BufferedReader br = new BufferedReader(new FileReader(fp))) {
			String line;
			while((line = br.readLine()) != null) {
				String[] tokens = line.trim().split(" "); //trim ws & tokenize
				processProblemLine(tokens);
			}
		} catch (IOException e) {
			System.err.format("IOException %s%n",e);
		}
	}

	public static void readOptTourLengths(String fp) {
		try(BufferedReader br = new BufferedReader(new FileReader(fp))) {
			String line;
			while((line = br.readLine()) != null) {
				String[] tokens = line.trim().split(" "); //trim ws & tokenize
				processOptTourLine(tokens);
			}
		} catch (IOException e) {
			System.err.format("IOException %s%n",e);
		}
		String problemName = problemFilePath.substring(
						problemFilePath.lastIndexOf('/') + 1,
						problemFilePath.lastIndexOf('.'));
		optimalTourCost = optTourLengths.get(problemName);
	}

	public static void readParams(String[] args) {
		if(args.length != 12 && args.length != 11) {
			printErrorAndExit();
		} 

		if(args[0].equals(acs)) {
			whichAlgorithm = 1;
			eps = Double.parseDouble(args[10]);
			qZero = Double.parseDouble(args[11]);
		} else if(args[0].equals(eas)) {
			whichAlgorithm = 0;
			e = Double.parseDouble(args[10]);
		} else {
			printErrorAndExit();
		}

		numberofAnts = Integer.parseInt(args[1]);
		numberofIterations = Integer.parseInt(args[2]);
		alpha = Double.parseDouble(args[3]);
		beta = Double.parseDouble(args[4]);
		rho = Double.parseDouble(args[5]);
		problemFilePath = args[6];
		stopCondition = Integer.parseInt(args[7]);
		secondsAllowed = Integer.parseInt(args[8]);
		errorAllowed = Double.parseDouble(args[9]);
	}

	public static String[] formatNodeInput(String[] s) {
		String t[] = new String[3];
		int counter = 0;
		for(int i = 0; i < s.length; i++) {
			if(!s[i].equals("")) {
				t[counter] = s[i];
				counter += 1;
			}
		}
		return t;
	}

	public static void outputSolution(ArrayList<Node> t, int c) {
		printSolution(t);
		System.out.println("Cost: " + c);
	}

	public static void printSolution(ArrayList<Node> s) {
		int lineCounter = 0;
		for(int i = 0; i < s.size(); i++) {
			System.out.print(s.get(i).getID() + "\t");
			if(lineCounter == 9) {
				System.out.println(); //ten variables per line
				lineCounter = 0;
			} else {
				lineCounter++;
			}
		}
		System.out.print("\n");
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
		System.out.println("Optimal Tour: " + optimalTourCost);
	}

	public static void printStringArray(String[] s) {
		for(int i = 0; i < s.length; i++) {
			System.out.print(s[i] + " ");
		}
		System.out.print("\n");
	}

	public static void printErrorAndExit() {
		System.out.println(incorrectParams);
		System.exit(1); //exit with error
	}
}