/*
ACO for TSP - Project 3
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

The code in this file contains the implementation of Ant Colony Optimization 
for the Traveling Salesman Problem, as part of Project 3. This file contains
implementations & logic for ACS and EAS ant systems.
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
	private static int secondsAllowed ;
	private static int optimalTourCost;
	private static double errorAllowed;
	
	//ACS parameters
	private static double eps;
	private static double tauZero; //calculated not provided
	private static double qZero;

	//EAS parameters
	private static double e;

	//Containers
	private static ArrayList<Node> nodes;
	private static ArrayList<Node> solutionTour;
	private static int solutionCost;
	private static double[][] pheromoneMatrix;
	private static Map<String, Integer> optTourLengths;
	private static String problemName;
	private static double[][] distanceMatrix;

	//Utility
	private static Random generator;

	//Time-keeping
	private static long timeStart;
	private static long timeFinish;

	//Visualizer
	Visualizer vis;

	public ACO() {
		nodes = new ArrayList<>();
		optTourLengths = new HashMap<String, Integer>();
		generator = new Random();
	}

	public void initializeACO(String[] args) {
		readParams(args);
		readProblem(problemFilePath);
		readOptTourLengths(Constants.optTourLengthsFilePath);
		distanceMatrix = initializeDistanceMatrix(nodes);
		initializeTauZero();
		printParams(whichAlgorithm);
	}

	public ArrayList<ArrayList<Node>> runACO() {

		timeStart = System.nanoTime();
		if(whichAlgorithm == 1) {
			solutionTour = acs();
		} else if(whichAlgorithm == 0) {
			solutionTour = eas();
		} else {
			Logger.printErrorAndExit(Constants.incorrectParams);
		}

		timeFinish = System.nanoTime();

		solutionCost = computeCost(solutionTour);
		outputSolution(solutionTour, solutionCost);
		// statsOutputSolution(solutionCost,((timeFinish - timeStart) / 1000000000.0));

		ArrayList<ArrayList<Node>> sol = new ArrayList<>();
		sol.add(nodes);
		sol.add(solutionTour);
		return sol;
	}

	//ACS function
	private static ArrayList<Node> acs() {
		//Create initial random solution
		ArrayList<Node> bestTour = initializeRandomSolution(nodes);
		int bestCost = computeCost(bestTour);
		//Initialize initial pheromone levels
		double initialPheromone = 1.0 / ((double)nodes.size() * bestCost);
		pheromoneMatrix = Utility.initializeMatrix(nodes.size(), initialPheromone);

		// //Initialize probability matrix used for storing already computed probabilities when determining the next city using
		// //the conventional city selection approach
		// double[][] probabilityMatrix;
		// //Initialize probability matrix used for storing already computed probabilities when determing the next city using
		// //the greedy heuristic approach 
		// double[][] acsProbabilityMatrix;

		ArrayList<Node> candidateTour;
		int candidateCost;

		int itCounter = 0;
		long st = System.nanoTime();

		while(evaluateStopCondition(stopCondition,itCounter, st, bestTour)) {
			itCounter += 1;
			
			for(int i = 0; i < numberofAnts; i++) {
				// System.out.println("TOUR START");
				candidateTour = constructSolutionACS();
				// System.out.println("TOUR END");
				candidateCost = computeCost(candidateTour);
				if(candidateCost < bestCost) {
					System.out.println("New Best: " + candidateCost + " at iteration " + itCounter);
					bestCost = candidateCost;
					bestTour = candidateTour;
				}

				localPheromoneUpdateACS(candidateTour);
			}
			// System.out.println("GLOBAL START");
			globalPheromoneUpdateACS();
			bestTourPheromoneUpdateACS(bestTour, bestCost);
			// System.out.println("GLOBAL END");
		}

		return bestTour;
	}

	private static ArrayList<Node> constructSolutionACS() {
		ArrayList<Integer> tour =  new ArrayList<>();
		int randomCity = generator.nextInt(nodes.size())+1;
		int nextCity;
		double r;

		tour.add(randomCity);

		while(tour.size() != nodes.size()) {
			 r = generator.nextDouble();
			 nextCity = (r <= qZero) ? pickNextCityACS(tour, tour.get(tour.size()-1)) : 
			 pickNextCityACSConventional(tour, tour.get(tour.size()-1));
			 tour.add(nextCity);
		}

		tour.add(randomCity);

		ArrayList<Node> rt = Utility.integerTourToNodeTour(tour,nodes);
		return rt;
	}

	private static int pickNextCityACS(ArrayList<Integer> tour, int lastCityID) {

		double currentMax = -Double.MAX_VALUE;
		int nextCityID = -1;
		double product;

		for(int i = 0; i < nodes.size(); i++) {
			if(!tour.contains(i+1)) {
				// if (acsProbabilityMatrix[lastCityID-1][i] == -1) {
					if (distanceMatrix[lastCityID-1][i] == 0.0) {
						return (i+1);
					}
					product = pheromoneMatrix[lastCityID-1][i] * 
					(Math.pow(1/distanceMatrix[lastCityID-1][i],beta));
				// }

				if (product > currentMax) {
					currentMax = product;
					nextCityID = i+1;
				}
			}
		}

		if(nextCityID == -1) {
			Logger.printErrorAndExit("Could not pick next city, ACS");
		}
		return nextCityID;
	}

	private static int pickNextCityACSConventional(ArrayList<Integer> tour, int lastCityID) {
		ArrayList<Double> nodeProbs = new ArrayList<>();
		ArrayList<Integer> nodeTracker = new ArrayList<>();
		double denomSum = 0.0;
		double product;

		for(int i = 0; i < nodes.size(); i++){
			if(!tour.contains(i+1)) {
				// if (probabilityMatrix[lastCityID-1][i] == -1) {
					if (distanceMatrix[lastCityID-1][i] == 0.0) {
						return (i+1);
					}
					product = (Math.pow(pheromoneMatrix[lastCityID-1][i], alpha)) * 
					(Math.pow(1/distanceMatrix[lastCityID-1][i], beta));
				// }
				nodeProbs.add(product);
				nodeTracker.add(i+1);
				denomSum += product;
			}
		}

		double cumulativeSum = 0;
		double random = generator.nextDouble();
		for (int i = 0; i < nodeProbs.size(); i++) {
			cumulativeSum += (nodeProbs.get(i)/denomSum);
			if (random < cumulativeSum) {
				if(tour.contains(nodeTracker.get(i))) {
					System.exit(1);
				}
				return nodeTracker.get(i);
			}
		}
		Logger.printErrorAndExit("Could not pick next city");
		return -1;
	}

	private static void localPheromoneUpdateACS(ArrayList<Node> cand) {
		int endIndex, startIndex;
		for(int i = 1; i < cand.size(); i++) {
			endIndex = cand.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = cand.get(i-1).getID() - 1;
			pheromoneMatrix[startIndex][endIndex]
						= ((1.0-eps) * pheromoneMatrix[startIndex][endIndex]) 
						+ (eps * tauZero);
		}
	}

	private static void globalPheromoneUpdateACS() {
		int matrixDimension = nodes.size();
		for (int i = 0; i < matrixDimension; i++) {
			for (int j = 0; j < matrixDimension; j++) {
				pheromoneMatrix[i][j] = ((1.0-rho)*pheromoneMatrix[i][j]);
			}
		}

	}

	private static void bestTourPheromoneUpdateACS(ArrayList<Node> bestTour, int bestCost) {
		int endIndex, startIndex;
		double updateValue = 1/bestCost;
		for(int i = 1; i < bestTour.size(); i++) {
			endIndex = bestTour.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = bestTour.get(i-1).getID() - 1;
			pheromoneMatrix[startIndex][endIndex] += (rho*updateValue);
		}
	}

	private static ArrayList<Node> eas() {
		ArrayList<Node> bestTour = initializeRandomSolution(nodes);
		int bestCost = computeCost(bestTour);
		double initialPheromone = 1.0 / ((double)nodes.size() * bestCost);
		pheromoneMatrix = Utility.initializeMatrix(nodes.size(), initialPheromone);

		double[][] legPheromoneUpdateMatrix;
		double[][] probabilityMatrix;
		 
		int itCounter = 0;
		long st = System.nanoTime();

		ArrayList<Node> candidateTour;
		int candidateCost;

		while(evaluateStopCondition(stopCondition,itCounter, st, bestTour)) {
			itCounter += 1;
			probabilityMatrix = Utility.initializeMatrix(nodes.size(),-1);
			legPheromoneUpdateMatrix = Utility.initializeMatrix(nodes.size(),0);

			for (int i = 0; i < numberofAnts; i++) {
				candidateTour = constructSolutionEAS(probabilityMatrix);
				candidateCost = computeCost(candidateTour);

				if (candidateCost < bestCost) {
					System.out.println("New Best: " + candidateCost + " at iteration " + itCounter);
					bestTour = candidateTour;
					bestCost = candidateCost;
				}
				legPheromoneUpdateEAS(candidateTour, candidateCost, legPheromoneUpdateMatrix);
			}
			pheromoneUpdateEAS(legPheromoneUpdateMatrix);
			bestTourPheromoneUpdate(bestTour, bestCost, e);
		}

		return bestTour;
	}

	private static ArrayList<Node> constructSolutionEAS(double[][] probabilityMatrix) {
		ArrayList<Integer> tour =  new ArrayList<>();
		int startCity = generator.nextInt(nodes.size()) + 1;
		tour.add(startCity); //Initial node
		int nextCity;

		while(tour.size() != nodes.size()) {
			nextCity = pickNextCity(tour, tour.get(tour.size()-1), probabilityMatrix);
			tour.add(nextCity);
		}
		tour.add(startCity);

		ArrayList<Node> rt = Utility.integerTourToNodeTour(tour,nodes);
		return rt;
	}

	private static void legPheromoneUpdateEAS(ArrayList<Node> candidateTour, int candidateCost,
										double[][] legPheromoneUpdateMatrix) {
		int endIndex, startIndex;
		double updateValue = 1/candidateCost;
		for(int i = 1; i < candidateTour.size(); i++) {
			endIndex = candidateTour.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = candidateTour.get(i-1).getID() - 1;
			legPheromoneUpdateMatrix[startIndex][endIndex] += updateValue;
		}
	}

	private static void pheromoneUpdateEAS(double[][] legPheromoneUpdateMatrix) {
		int matrixDimension = nodes.size();

		for (int i = 0; i < matrixDimension; i++) {
			for (int j = 0; j < matrixDimension; j++) {
				pheromoneMatrix[i][j] = ((1-rho)
										* pheromoneMatrix[i][j]) 
										+ legPheromoneUpdateMatrix[i][j];
			}
		}
	}

	private static void bestTourPheromoneUpdate(ArrayList<Node> bestTour, int bestCost, double e) {
		int endIndex, startIndex;
		double updateValue = 1/bestCost;
		for(int i = 1; i < bestTour.size(); i++) {
			endIndex = bestTour.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = bestTour.get(i-1).getID() - 1;
			pheromoneMatrix[startIndex][endIndex] += (e*updateValue);
		}
	}

	private static int pickNextCity(ArrayList<Integer> tour, int lastCityID, double[][] probabilityMatrix) {
		ArrayList<Double> nodeProbs = new ArrayList<>();
		ArrayList<Integer> nodeTracker = new ArrayList<>();
		double denomSum = 0.0;

		for(int i = 0; i < nodes.size(); i++){
			if(!tour.contains(i+1)) {
				if (probabilityMatrix[lastCityID-1][i] == -1) {
					if (distanceMatrix[lastCityID-1][i] == 0.0) {
						return (i+1);
					}
					probabilityMatrix[lastCityID-1][i] = (Math.pow(pheromoneMatrix[lastCityID-1][i], alpha)) * 
					(Math.pow(1/distanceMatrix[lastCityID-1][i], beta));
				}
				nodeProbs.add(probabilityMatrix[lastCityID-1][i]);
				nodeTracker.add(i+1);
				denomSum += probabilityMatrix[lastCityID-1][i];
			}
		}

		double cumulativeSum = 0;
		double random = generator.nextDouble();
		for (int i = 0; i < nodeProbs.size(); i++) {
			cumulativeSum += (nodeProbs.get(i)/denomSum);
			if (random < cumulativeSum) {
				if(tour.contains(nodeTracker.get(i))) {
					System.exit(1);
				}
				return nodeTracker.get(i);
			}
		}
		Logger.printErrorAndExit("Could not pick next city");
		return -1;
	}

	// 0 - terminate after max iterations reached
	// 1 - terminate if found tour that is no more than a specified percentage 
	//	over the optimal (0.0 would mean you will not settle for anything less than the optimal)
	// 2 - both
	// 3 - terminate after secondsAllowed exceeds
	// 4 - all three
	private static Boolean evaluateStopCondition(int stopCondition, int currIt, 
										long startTime, ArrayList<Node> tour) {

		switch(stopCondition) {
			case 0: 
				System.out.print("Iteration: " + currIt + "\r");
				return (currIt < numberofIterations);
			case 1:	
				System.out.println("Error: " + ((computeCost(tour) / (double)(optimalTourCost))) + "\r");
				return (((computeCost(tour) / (double)(optimalTourCost))) > (errorAllowed + 1));
			case 2:
				System.out.print("Seconds Passed: " + ((System.nanoTime() - startTime) / 1000000000.0) + "\r");
				return (((System.nanoTime() - startTime) / 1000000000.0) < secondsAllowed);
			case 3:
				System.out.print("Error: " + ((computeCost(tour) / (double)(optimalTourCost))) + "\t"
					+ "Seconds Passed: " + ((System.nanoTime() - startTime) / 1000000000.0) + "\r");
				return ((((computeCost(tour) / (double)(optimalTourCost))) > (errorAllowed + 1))
					|| (currIt < numberofIterations));
			default:
				Logger.printErrorAndExit("Unknown stop condition");
		}
		return false;
	}

	private static ArrayList<Node> initializeRandomSolution(ArrayList<Node> s) {
		ArrayList<Node> randSol = s;
		long seed = System.nanoTime();
		Collections.shuffle(randSol, new Random(seed));
		return randSol;
	}

	private static void initializeTauZero() {
		ArrayList<Integer> nearestNeighborTour = new ArrayList<>();
		int randomCity = generator.nextInt(nodes.size())+1;
		nearestNeighborTour.add(randomCity);

		while(nearestNeighborTour.size() != nodes.size()) {
			int closestNeighbor = Utility.getClosestCityTo(
									nearestNeighborTour.get(
											nearestNeighborTour.size()-1),
									nearestNeighborTour,nodes);
			nearestNeighborTour.add(closestNeighbor);
		}

		ArrayList<Node> rt = Utility.integerTourToNodeTour(nearestNeighborTour,nodes);
		rt.add(Utility.getCity(nearestNeighborTour.get(0),nodes)); //add first city as last one to complete tour

		tauZero = (1.0/(numberofAnts * computeCost(rt)));
	}

	private static double[][] initializeDistanceMatrix(ArrayList<Node> nodes) {
		double[][] temp = new double[nodes.size()][nodes.size()];

		for(int r = 0; r < nodes.size(); r++) {
			for(int c = 0; c < nodes.size(); c++) {
				temp[r][c] = temp[c][r] = Utility.euclideanDistance2D(Utility.getCity(r+1,nodes), Utility.getCity(c+1,nodes));
			}
		}
		return temp;
	}

	//compute cost of the given tour by adding up the distance between successive nodes
	public static int computeCost(ArrayList<Node> t) {
		int distance = 0;
		for(int i = 1; i < t.size(); i++) {
			distance += distanceMatrix[t.get(i).getID()-1][t.get(i-1).getID()-1];
		}
		return distance;
	}

	private static void processProblemLine(String[] tokens) {
		try {
			Integer.parseInt(tokens[0]); //make sure its a node line
			tokens = Utility.formatNodeInput(tokens);

			Node temp = new Node(Integer.parseInt(tokens[0]),
							Double.parseDouble(tokens[1]),
							Double.parseDouble(tokens[2]));
			nodes.add(temp);
		} catch (NumberFormatException e) {
			// System.err.format("NumberFormatException %s%n",e);
		}
	}

	private static void processOptTourLine(String[] tokens) {
		try {
			optTourLengths.put(tokens[0],Integer.parseInt(tokens[1]));
		} catch (NumberFormatException e) {
			// System.err.format("NumberFormatException %s%n",e);
		}
	}

	private static void readProblem(String fp) {
		try(BufferedReader br = new BufferedReader(new FileReader(fp))) {
			String line;
			while((line = br.readLine()) != null) {
				String[] tokens = line.trim().split(" "); //trim ws & tokenize
				processProblemLine(tokens);
			}
		} catch (IOException e) {
			System.err.format("IOException %s%n",e);
			Logger.printErrorAndExit("Problem File Not Found. Exiting");
		}
	}

	private static void readOptTourLengths(String fp) {
		try(BufferedReader br = new BufferedReader(new FileReader(fp))) {
			String line;
			while((line = br.readLine()) != null) {
				String[] tokens = line.trim().split(" "); //trim ws & tokenize
				processOptTourLine(tokens);
			}
		} catch (IOException e) {
			System.err.format("IOException %s%n",e);
		}
		problemName = problemFilePath.substring(
						problemFilePath.lastIndexOf('/') + 1,
						problemFilePath.lastIndexOf('.'));
		optimalTourCost = optTourLengths.get(problemName);
	}

	private static void readParams(String[] args) {
		if(args.length != 13 && args.length != 12) {
			Logger.printErrorAndExit(Constants.incorrectParams);
		} 

		if(args[0].equals(Constants.acs)) {
			whichAlgorithm = 1;
			eps = Double.parseDouble(args[10]);
			qZero = Double.parseDouble(args[11]);
		} else if(args[0].equals(Constants.eas)) {
			whichAlgorithm = 0;
			e = Double.parseDouble(args[10]);
		} else {
			Logger.printErrorAndExit("Unknown algorithm");
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

	private static void outputSolution(ArrayList<Node> t, int c) {
		System.out.print("\n");
		Logger.printNodeArrayList(t);
		System.out.println("Cost: " + c);
		System.out.println("Error: " + ((1.0*c)/optimalTourCost));
	}

	private static void statsOutputSolution(int c, double secElapsed) {
		System.out.println(problemName + " " + optimalTourCost + " " + c + " " 
								+ ((1.0*c)/optimalTourCost) + " " + secElapsed);
	}

	private static void printParams(int alg) {
		System.out.println("Algorithm: " + ((alg == 1) ? Constants.acs : Constants.eas)
			+ "\n# of Ants: " + numberofAnts
			+ "\n# of Iterations: " + numberofIterations
			+ "\n# of Cities: " + nodes.size()
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
}