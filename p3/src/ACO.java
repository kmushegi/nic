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

	//Utility
	private static Random generator;

	public ACO() {
		nodes = new ArrayList<>();
		optTourLengths = new HashMap<String, Integer>();
		generator = new Random();
	}

	public void initializeACO(String[] args) {
		readParams(args);
		readProblem(problemFilePath);
		readOptTourLengths(Constants.optTourLengthsFilePath);
		initializeTauZero();
		printParams(whichAlgorithm);
	}

	public void runACO() {
		if(whichAlgorithm == 1) {
			solutionTour = acs(numberofAnts, numberofIterations, alpha, beta, 
								rho, eps, qZero);
		} else if(whichAlgorithm == 0) {
			solutionTour = eas(numberofAnts, numberofIterations, alpha, beta, 
								rho, e);
		} else {
			Logger.printErrorAndExit(Constants.incorrectParams);
		}

		solutionCost = Utility.computeCost(solutionTour);
		outputSolution(solutionTour, solutionCost);
	}

	private static ArrayList<Node> acs(int ants, int its, double a, double b,
									double rho, double eps, double q) {
		ArrayList<Node> bestTour = initializeRandomSolution(nodes);
		int bestCost = Utility.computeCost(bestTour);
		double initialPheromone = 1.0 / ((double)nodes.size() * bestCost);
		pheromoneMatrix = initializePheromoneMatrix(nodes.size(), initialPheromone);

		int itCounter = 0;
		long st = System.nanoTime();

		while(evaluateStopCondition(stopCondition,itCounter, st, bestTour)) {
			itCounter += 1;
			
			for(int i = 0; i < ants; i++) {
				ArrayList<Node> candidateTour;
				int candidateCost;

				candidateTour = constructSolutionACS(pheromoneMatrix, b, q);
				candidateCost = Utility.computeCost(candidateTour);
				if(candidateCost < bestCost) {
					bestCost = candidateCost;
					bestTour = new ArrayList<>(candidateTour);
				}
				localPheromoneUpdateACS(candidateTour, eps);
			}
			globalPheromoneUpdateACS(rho);
			bestTourPheromoneUpdateACS(bestTour, rho);
		}

		return bestTour;
	}

	private static ArrayList<Node> constructSolutionACS(double[][] pm, double b, double q) {
		ArrayList<Integer> tour =  new ArrayList<>();
		int randomCity = generator.nextInt(nodes.size())+1;
		tour.add(randomCity);

		while(tour.size() != nodes.size()) {
			int nextCity = pickNextCityACS(tour,q,tour.get(tour.size()-1));
			tour.add(nextCity);
		}

		tour.add(randomCity);

		ArrayList<Node> rt = Utility.integerTourToNodeTour(tour,nodes);
		return rt;
	}

	private static int pickNextCityACS(ArrayList<Integer> tour, double q, int lastCityID) {
		double r = generator.nextDouble();

		double currentMax = -Double.MAX_VALUE;
		int nextCityID = -1;

		if(r <= q) {
			for(int i = 0; i < nodes.size(); i++) {
				if(!tour.contains(i+1)) {
					double invDistance = 1/(Utility.euclideanDistance2D(
											Utility.getCity(lastCityID,nodes), 
											Utility.getCity(i+1,nodes)));
					double pLevel = pheromoneMatrix[lastCityID-1][i];

					double temp = pLevel * (Math.pow(invDistance,beta));
					if(temp > currentMax) {
						currentMax = temp;
						nextCityID = (i+1);
					}
				}
			}
		} else {
			nextCityID = pickNextCityEAS(tour, beta, alpha, tour.get(tour.size()-1));
		}
		if(nextCityID == -1) {
			Logger.printErrorAndExit("Could not pick next city, ACS");
		}
		return nextCityID;
	}

	private static void localPheromoneUpdateACS(ArrayList<Node> cand, double eps) {
		int endIndex, startIndex;
		for(int i = 1; i < cand.size(); i++) {
			endIndex = cand.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = cand.get(i-1).getID() - 1;
			pheromoneMatrix[startIndex][endIndex]
						= ((1.0-eps) * pheromoneMatrix[startIndex][endIndex]) 
						+ (eps * tauZero);
		}
	}

	private static void globalPheromoneUpdateACS(double rho) {
		int matrixDimension = nodes.size();
		for (int i = 0; i < matrixDimension; i++) {
			for (int j = 0; j < matrixDimension; j++) {
				pheromoneMatrix[i][j] = ((1.0-rho)*pheromoneMatrix[i][j]);
				pheromoneMatrix[j][i] = ((1.0-rho)*pheromoneMatrix[j][i]);
			}
		}

	}

	private static void bestTourPheromoneUpdateACS(ArrayList<Node> bestTour, double rho) {
		int endIndex, startIndex;
		double updateValue;
		for(int i = 1; i < bestTour.size(); i++) {
			updateValue = 1/(Utility.euclideanDistance2D(bestTour.get(i), bestTour.get(i-1)));
			endIndex = bestTour.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = bestTour.get(i-1).getID() - 1;
			pheromoneMatrix[startIndex][endIndex] += (rho*updateValue);
		}
	}

	private static ArrayList<Node> eas(int ants, int its, double a, double b,
									double rho, double e) {

		ArrayList<Node> bestTour = initializeRandomSolution(nodes);
		int bestCost = Utility.computeCost(bestTour);
		double initialPheromone = 1.0 / ((double)nodes.size() * bestCost);
		pheromoneMatrix = initializePheromoneMatrix(nodes.size(), initialPheromone);

		double[][] legPheromoneUpdateMatrix = initializeLegPheromoneUpdateMatrix(nodes.size());
		 
		int itCounter = 0;
		long st = System.nanoTime();

		ArrayList<Node> candidateTour;
		int candidateCost;

		while(evaluateStopCondition(stopCondition,itCounter, st, bestTour)) {
			itCounter += 1;

			for (int i = 0; i < ants; i++) {
				candidateTour = constructSolutionEAS(b, a, rho);
				candidateCost = Utility.computeCost(candidateTour);

				if (candidateCost < bestCost) {
					System.out.println("New Best: " + candidateCost);

					bestTour = candidateTour;
					bestCost = candidateCost;
				}
				legPheromoneUpdateEAS(candidateTour, legPheromoneUpdateMatrix);
			}
			pheromoneUpdateEAS(rho, legPheromoneUpdateMatrix);
			bestTourPheromoneUpdate(bestTour, e);
		}

		return bestTour;
	}

	private static ArrayList<Node> constructSolutionEAS(double b, double a, double rho) {
		ArrayList<Integer> tour =  new ArrayList<>();
		int startCity = generator.nextInt(nodes.size()) + 1;
		tour.add(startCity); //Initial node
		int nextCity;

		while(tour.size() != nodes.size()) {
			nextCity = pickNextCityEAS(tour, b, a, tour.get(tour.size()-1));
			tour.add(nextCity);
		}
		tour.add(startCity);

		ArrayList<Node> rt = Utility.integerTourToNodeTour(tour,nodes);
		return rt;
	}

	private static int pickNextCityEAS(ArrayList<Integer> tour, double b, double a, int lastCityID) {
		ArrayList<Double> nodeProbs = new ArrayList<>();
		ArrayList<Integer> nodeTracker = new ArrayList<>();
		double denomSum = 0.0;

		// System.out.println("Last City ID: " + lastCityID);

		for(int i = 0; i < nodes.size(); i++){
			if(!tour.contains(i+1)) {
				// System.out.println("Curr City I: " + (i+1));
				double distance = Utility.euclideanDistance2D(
										Utility.getCity(lastCityID,nodes), 
										Utility.getCity(i+1,nodes));

				if (distance == 0.0) {
					return i+1;
				}

				double invDistance = 1/(distance);
				// System.out.println("Dist: " + invDistance);
				double pLevel = pheromoneMatrix[lastCityID-1][i];
				double nodeNum = (Math.pow(pLevel, a)) * (Math.pow(invDistance, b));
				nodeProbs.add(nodeNum);
				nodeTracker.add(i+1);
				denomSum += nodeNum;
			}
		}

		for (int i = 0; i < nodeProbs.size(); i++){
			nodeProbs.set(i, (nodeProbs.get(i)/denomSum));
		}

		// double sum = 0;
		// for (int i = 0 ; i < nodeProbs.size(); i++) {
		// 	sum += nodeProbs.get(i);
		// }

		// System.out.println("Total Prob: " + sum);

		double cumulativeSum = 0;
		double random = generator.nextDouble();
		for (int i = 0; i < nodeProbs.size(); i++) {
			cumulativeSum += nodeProbs.get(i);
			if (random < cumulativeSum) {
				if(tour.contains(nodeTracker.get(i))) {
					System.exit(1);
				}
				return nodeTracker.get(i);
			}
		}
		Logger.printErrorAndExit("Could not pick next city, EAS");
		return -1;
	}

	private static void legPheromoneUpdateEAS(ArrayList<Node> candidateTour, 
										double[][] legPheromoneUpdateMatrix) {
		int endIndex, startIndex;
		double updateValue;

		for(int i = 1; i < candidateTour.size(); i++) {
			updateValue = 1/(Utility.euclideanDistance2D(candidateTour.get(i), candidateTour.get(i-1)));
			endIndex = candidateTour.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = candidateTour.get(i-1).getID() - 1;
			legPheromoneUpdateMatrix[startIndex][endIndex] += updateValue;
		}
	}

	private static void pheromoneUpdateEAS(double rho, double[][] legPheromoneUpdateMatrix) {
		int matrixDimension = nodes.size();

		for (int i = 0; i < matrixDimension; i++) {
			for (int j = 0; j < matrixDimension; j++) {
				pheromoneMatrix[i][j] = ((1-rho)
										* pheromoneMatrix[i][j]) 
										+ legPheromoneUpdateMatrix[i][j];
				pheromoneMatrix[j][i] = ((1-rho)
										* pheromoneMatrix[j][i]) 
										+ legPheromoneUpdateMatrix[j][i];
			}
		}
	}

	private static void bestTourPheromoneUpdate(ArrayList<Node> bestTour, double e) {
		int endIndex, startIndex;
		double updateValue;
		for(int i = 1; i < bestTour.size(); i++) {
			updateValue = 1/(Utility.euclideanDistance2D(bestTour.get(i), bestTour.get(i-1)));
			endIndex = bestTour.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startIndex = bestTour.get(i-1).getID() - 1;
			pheromoneMatrix[startIndex][endIndex] += (e*updateValue);
		}
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
				System.out.print("Error: " + ((Utility.computeCost(tour) / optimalTourCost) - 1) + "\r");
				return (((Utility.computeCost(tour) / optimalTourCost) - 1) > errorAllowed);
			// case 2:
			// 	return ((currIt < numberofIterations) 
			// 		&& (((Utility.computeCost(tour) / optimalTourCost) - 1) > errorAllowed));
			case 2:
				System.out.print("Seconds Passed: " + ((System.nanoTime() - startTime) / 1000000000.0) + "\r");
				return (((System.nanoTime() - startTime) / 1000000000.0) < secondsAllowed);
			case 3:
				System.out.print("Error: " + ((Utility.computeCost(tour) / optimalTourCost) - 1) + "\t"
					+ "Seconds Passed: " + ((System.nanoTime() - startTime) / 1000000000.0) + "\r");
				return ((((Utility.computeCost(tour) / optimalTourCost) - 1) > errorAllowed)
					&& ((System.nanoTime() - startTime) / 1000000000.0) < secondsAllowed);
			default:
				Logger.printErrorAndExit("Unknown stop condition");
		}
		return false;
	}

	private static double[][] initializePheromoneMatrix(int n, double initPh) {
		double[][] temp = new double[n][n];

		for(int r = 0; r < n; r++) {
			for(int c = 0; c < n; c++) {
				temp[r][c] = initPh;
			}
		}
		return temp;
	}

	private static double[][] initializeLegPheromoneUpdateMatrix(int n) {
		double[][] temp = new double[n][n];

		for(int r = 0; r < n; r++) {
			for(int c = 0; c < n; c++) {
				temp[r][c] = 0;
			}
		}
		return temp;
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
			int closestNeighbor = Utility.getClosestCityTo(nearestNeighborTour.get(
										nearestNeighborTour.size()-1),nearestNeighborTour,nodes);
			nearestNeighborTour.add(closestNeighbor);
		}

		ArrayList<Node> rt = Utility.integerTourToNodeTour(nearestNeighborTour,nodes);
		rt.add(Utility.getCity(nearestNeighborTour.get(0),nodes)); //add first city as last one to complete tour

		tauZero = (1.0/(numberofAnts * Utility.computeCost(rt)));
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
		String problemName = problemFilePath.substring(
						problemFilePath.lastIndexOf('/') + 1,
						problemFilePath.lastIndexOf('.'));
		optimalTourCost = optTourLengths.get(problemName);
	}

	private static void readParams(String[] args) {
		if(args.length != 12 && args.length != 11) {
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