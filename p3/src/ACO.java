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

		ArrayList<Node> candidateTour;
		int candidateCost;

		int itCounter = 0;
		long st = System.nanoTime();

		//While the stop condition has not been reached, continue to generate solutions
		while(evaluateStopCondition(stopCondition,itCounter, st, bestTour)) {
			itCounter += 1;
			
			//For all antss
			for(int i = 0; i < numberofAnts; i++) {
				//Construct a candidate solution
				candidateTour = constructSolutionACS();
				candidateCost = computeCost(candidateTour);

				//If the candidate's tour cost is less that the current best, it becomes the new best tour
				if(candidateCost < bestCost) {
					System.out.println("New Best: " + candidateCost + " at iteration " + itCounter);
					bestCost = candidateCost;
					bestTour = candidateTour;
				}

				//Perform the local pheromone update on the most recent tour
				localPheromoneUpdateACS(candidateTour);
			}

			//Perform the global pheromone update on all legs between cities
			globalPheromoneUpdateACS();
			//Perform the pheromone addition to the best current tour
			bestTourPheromoneUpdateACS(bestTour, bestCost);
		}

		return bestTour;
	}

	private static ArrayList<Node> constructSolutionACS() {
		ArrayList<Integer> tour =  new ArrayList<>();
		//Generate a random starting city ID
		int randomCity = generator.nextInt(nodes.size())+1;
		int nextCity;
		double r;

		//Add the random start ID to the current candidate tour
		tour.add(randomCity);

		//While not all cities have been visited, continue to populate tour with cities
		while(tour.size() != nodes.size()) {
			 r = generator.nextDouble();
			 //With probability qzero, perform the greedy heuristic, otherwise, conventionally select the next city
			 nextCity = (r <= qZero) ? pickNextCityACS(tour, tour.get(tour.size()-1)) : 
			 pickNextCityACSConventional(tour, tour.get(tour.size()-1));
			 //Add the new city to the tour
			 tour.add(nextCity);
		}

		//Add the starting city to the tour again as for TSP, we need to return to the starting city
		tour.add(randomCity);

		ArrayList<Node> rt = Utility.integerTourToNodeTour(tour,nodes);
		return rt;
	}

	private static int pickNextCityACS(ArrayList<Integer> tour, int lastCityID) {

		double currentMax = -Double.MAX_VALUE;
		int nextCityID = -1;
		double product;

		for(int i = 0; i < nodes.size(); i++) {
			//Pick a city that is not already in the tour
			if(!tour.contains(i+1)) {
				//If the distance to the city is zero, that is the same city (accounting for problems with two instances of the same city)
				if (distanceMatrix[lastCityID-1][i] == 0.0) {
					return (i+1);
				}
				//Calculate the product of pheromone levels * inverse distance^beta for the current city and new city pair
				product = pheromoneMatrix[lastCityID-1][i] * 
				(Math.pow(1/distanceMatrix[lastCityID-1][i],beta));

				//If product is less than the current max, product becomes the new max
				if (product > currentMax) {
					currentMax = product;
					nextCityID = i+1;
				}
			}
		}

		if(nextCityID == -1) {
			Logger.printErrorAndExit("Could not pick next city, ACS");
		}
		//Return the city that maximizes pheromone levels * inverse distance^beta
		return nextCityID;
	}

	private static int pickNextCityACSConventional(ArrayList<Integer> tour, int lastCityID) {
		ArrayList<Double> nodeProbs = new ArrayList<>();
		ArrayList<Integer> nodeTracker = new ArrayList<>();
		double denomSum = 0.0;
		double product;

		for(int i = 0; i < nodes.size(); i++){
			//Pick a city that hasn't been visited yet
			if(!tour.contains(i+1)) {
				if (distanceMatrix[lastCityID-1][i] == 0.0) {
					return (i+1);
				}
				//Compute the product of pheromone levels^alpha * inverse distance^beta
				product = (Math.pow(pheromoneMatrix[lastCityID-1][i], alpha)) * 
				(Math.pow(1/distanceMatrix[lastCityID-1][i], beta));
				//Add product to probability array
				nodeProbs.add(product);
				//Add city instance to city array
				nodeTracker.add(i+1);
				//Add product to the sum that will normalize probabilities
				denomSum += product;
			}
		}

		double cumulativeSum = 0;
		double random = generator.nextDouble();

		//Probabilistically pick the next city
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
			//Perform the pheromone update on all the legs of the most recent tour
			endIndex = cand.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = cand.get(i-1).getID() - 1;
			pheromoneMatrix[startIndex][endIndex]
						= ((1.0-eps) * pheromoneMatrix[startIndex][endIndex]) 
						+ (eps * tauZero);
		}
	}

	private static void globalPheromoneUpdateACS() {
		int matrixDimension = nodes.size();
		//Perform the global pheromone update on all legs between all cities in the graph
		for (int i = 0; i < matrixDimension; i++) {
			for (int j = 0; j < matrixDimension; j++) {
				pheromoneMatrix[i][j] = ((1.0-rho)*pheromoneMatrix[i][j]);
			}
		}
	}

	private static void bestTourPheromoneUpdateACS(ArrayList<Node> bestTour, int bestCost) {
		int endIndex, startIndex;
		double updateValue = 1/bestCost;
		//Add additional pheromone on the legs of the current best tour so as to favor the current best tour
		for(int i = 1; i < bestTour.size(); i++) {
			endIndex = bestTour.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = bestTour.get(i-1).getID() - 1;
			pheromoneMatrix[startIndex][endIndex] += (rho*updateValue);
		}
	}

	//EAS function
	private static ArrayList<Node> eas() {
		ArrayList<Node> bestTour = initializeRandomSolution(nodes);
		int bestCost = computeCost(bestTour);
		double initialPheromone = 1.0 / ((double)nodes.size() * bestCost);
		pheromoneMatrix = Utility.initializeMatrix(nodes.size(), initialPheromone);

		//Initialize matrix to store the pheromone deposited on each leg by the total number of ants traversing that leg
		double[][] legPheromoneUpdateMatrix;
		//Initialize the matricies used to store already calculated probability values used when picking the next city
		double[][] probabilityMatrix;
		 
		int itCounter = 0;
		long st = System.nanoTime();

		ArrayList<Node> candidateTour;
		int candidateCost;

		while(evaluateStopCondition(stopCondition,itCounter, st, bestTour)) {
			itCounter += 1;
			//Reset matrix as the ants that traverse each leg change from iteration to iteration
			legPheromoneUpdateMatrix = Utility.initializeMatrix(nodes.size(),0);
			//Reset probability matrix, as the pheromones on which probabilities are based change from iteration to iteration
			//after the pheromone updates
			probabilityMatrix = Utility.initializeMatrix(nodes.size(),-1);

			for (int i = 0; i < numberofAnts; i++) {
				candidateTour = constructSolutionEAS(probabilityMatrix);
				candidateCost = computeCost(candidateTour);

				if (candidateCost < bestCost) {
					System.out.println("New Best: " + candidateCost + " at iteration " + itCounter);
					bestTour = candidateTour;
					bestCost = candidateCost;
				}
				//Add pheromone to the matrix storing the pheromone deposited on each leg by the total number of ants traversing that leg
				legPheromoneUpdateEAS(candidateTour, candidateCost, legPheromoneUpdateMatrix);
			}
			//Update the pheromone on each leg of the graph
			pheromoneUpdateEAS(legPheromoneUpdateMatrix);
			//Add additional pheromone on the best so far tour
			bestTourPheromoneUpdate(bestTour, bestCost, e);
		}

		return bestTour;
	}

	private static ArrayList<Node> constructSolutionEAS(double[][] probabilityMatrix) {
		ArrayList<Integer> tour =  new ArrayList<>();
		int startCity = generator.nextInt(nodes.size()) + 1;
		tour.add(startCity);
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
		//Add pheromone to the matrix storing the pheromone deposited by all ants traversing specific leg
		for(int i = 1; i < candidateTour.size(); i++) {
			endIndex = candidateTour.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = candidateTour.get(i-1).getID() - 1;
			legPheromoneUpdateMatrix[startIndex][endIndex] += updateValue;
		}
	}

	private static void pheromoneUpdateEAS(double[][] legPheromoneUpdateMatrix) {
		int matrixDimension = nodes.size();

		//Perform pheromone evaporation, and add the pheromone deposited by all ants traversing specific leg
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
		//Add additional pheromone on the legs of the current best tour so as to favor the current best tour
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
	// 2 - terminate after secondsAllowed exceeds
	// 3 - terminate after 1 or 0 is satisfied
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

	//Randomly shuffle cities to intialize a random solution
	private static ArrayList<Node> initializeRandomSolution(ArrayList<Node> s) {
		ArrayList<Node> randSol = s;
		long seed = System.nanoTime();
		Collections.shuffle(randSol, new Random(seed));
		return randSol;
	}

	//Calulate nearest nieghbor for current problem
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

	//Pre-calculate the distances between cities in the current problem
	private static double[][] initializeDistanceMatrix(ArrayList<Node> nodes) {
		double[][] temp = new double[nodes.size()][nodes.size()];

		for(int r = 0; r < nodes.size(); r++) {
			for(int c = 0; c < nodes.size(); c++) {
				temp[r][c] = temp[c][r] = Utility.euclideanDistance2D(
											Utility.getCity(r+1,nodes), 
											Utility.getCity(c+1,nodes));
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