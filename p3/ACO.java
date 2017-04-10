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

	//MARCUS: Stop conditions. I don't think we can do conditions 0 and x together, because then it would basically just be max iterations

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
		// initializeTauZero();

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
			System.out.println("Iteration: " + itCounter);
			//store solutions here?
			for(int i = 0; i < ants; i++) {
				ArrayList<Node> candidateTour;
				int candidateCost;

				candidateTour = constructSolutionACS(pheromoneMatrix, b, q);
				candidateCost = computeCost(candidateTour);
				if(candidateCost < bestCost) {
					bestCost = candidateCost;
					bestTour = new ArrayList<>(candidateTour);
				}
				localPheromoneUpdate(candidateTour, t, initialPheromone);
			}
			globalPheromoneUpdate(bestTour, a);
		}

		return bestTour;
	}

	public static ArrayList<Node> constructSolutionACS(double[][] pm, double b, double q) {
		ArrayList<Integer> tour =  new ArrayList<>();
		int randomCity = generator.nextInt(nodes.size())+1;
		// System.out.println("Random City: " + randomCity);
		tour.add(randomCity);

		while(tour.size() != nodes.size()) {
			// if(tour.get(tour.size()-1) == 0) {
			// 	System.out.println("ALRT");
			// 	System.exit(1);
			// }

			ArrayList<Map<String,Double>> ch = generateChoices(tour.get(tour.size()-1), tour, b, q, 1.0);

			// System.out.println("Number of Choices: " + ch.size());

			// for (int i = 0; i < ch.size(); i++) {
			// 	System.out.println(ch.get(i).get("city"));
			// }

			int nextCity = pickNextCityACS(ch);
			// System.out.println(nextCity);
			// if(nextCity == 0) {
			// 	System.out.println("ALRT3");
			// 	System.exit(1);
			// }
			tour.add(nextCity);
		}

		tour.add(randomCity);

		//convert Integer tour to Node tour.
		ArrayList<Node> rt = new ArrayList<>();
		for(int i = 0; i < tour.size(); i++) {
			rt.add(getCity(tour.get(i)));
		}
		return rt;
	}

	public static ArrayList<Map<String,Double>> generateChoices(int lastCity, 
					ArrayList<Integer> tour, double b, double q, double hist) {
		ArrayList<Map<String, Double>> ch = new ArrayList<>();

		for(int i = 0; i < nodes.size(); i++) {
			if(!tour.contains((i+1))) {
				Map<String, Double> temp = new HashMap<String, Double>();
				temp.put("city",(double)(i+1));
				temp.put("hist",Math.pow(pheromoneMatrix[lastCity-1][i],hist));
				temp.put("dist",euclideanDistance2D(getCity(lastCity),getCity(i+1)));
				temp.put("heur",Math.pow((1.0/temp.get("dist")),q));
				temp.put("prob",(temp.get("hist") * temp.get("heur")));
				ch.add(temp);
			}
		}
		return ch;
	}

	public static int pickNextCityACS(ArrayList<Map<String,Double>> ch) {
		double sum = 0.0;
		for(int i = 0; i < ch.size(); i++) {
			sum += ch.get(i).get("prob");
			// System.out.println("adding: " + ch.get(i).get("prob"));
			// System.out.println("Sum: " + sum);
		}
		// System.exit(1);
		if(sum == 0.0) {
			// System.out.println("Next City 1: " + ch.get((generator.nextInt(nodes.size())+1)).get("city"));
			return ch.get((generator.nextInt(nodes.size())+1)).get("city").intValue();
		}

		double r = generator.nextDouble();
		for(int i = 0; i < ch.size(); i++) {
			r -= (ch.get(i).get("prob")/sum);
			if(r <= 0.0) {
				// System.out.println("Next City 2: " + ch.get(i).get("city"));
				return ch.get(i).get("city").intValue();
			}
		}
		// System.out.println("Next City 3: " + ch.get(ch.size()-1).get("city").intValue());
		return ch.get(ch.size()-1).get("city").intValue();
	}

	public static void localPheromoneUpdate(ArrayList<Node> cand, double t, double initPhero) {

		// System.out.println("Nodes Size: " + nodes.size());
		// System.out.println("Cands Size: " + cand.size());

		// for (int i = 0; i < cand.size(); i++) {
		// 	System.out.println(cand.get(i).getID());
		// }

		// System.exit(0);

		int endID, startID;
		for(int i = 1; i < cand.size(); i++) {
			endID = cand.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startID = cand.get(i-1).getID() - 1;

			double value = ((1.0-t)*pheromoneMatrix[startID][endID])+(t*initPhero);
			pheromoneMatrix[startID][endID] = value;
			pheromoneMatrix[endID][startID] = value;
		}
	}

	public static void globalPheromoneUpdate(ArrayList<Node> best, double a) {
		int endID, startID;
		for(int i = 1; i < best.size(); i++) {
			endID = best.get(i).getID() - 1; //Minus 1 as matrix is zero indexed
			startID = best.get(i-1).getID() - 1;

			double value = ((1.0-a)*pheromoneMatrix[startID][endID])+(a*(1.0/computeCost(best)));
			pheromoneMatrix[startID][endID] = value;
			pheromoneMatrix[endID][startID] = value;
		}
	}

	//WHAT IS E?
	public static ArrayList<Node> eas(int ants, int its, double a, double b,
									double rho, double e) {
		ArrayList<Node> bestTour = initializeRandomSolution(nodes);
		int bestCost = computeCost(bestTour);
		double initialPheromone = 1.0 / ((double)nodes.size() * bestCost);
		pheromoneMatrix = initializePheromoneMatrix(nodes.size(), initialPheromone);

		int itCounter = 0;
		long st = System.nanoTime();

		while(evaluateStopCondition(stopCondition,itCounter, st, bestTour)) {
			itCounter += 1;

			for (int i = 0; i < ants; i++) {
				ArrayList<Node> candidateTour = constructSolutionEAS(pheromoneMatrix, b , a, rho);
				int candidateCost = computeCost(candidateTour);

				if (candidateCost < bestCost) {
					bestTour = candidateTour;
					bestCost = candidateCost;
				}
			}

			//globalPheromoneUpdateEAS(candidateTour , bestTour , pheromoneMatrix);

		}

		return bestTour;
	}

	public static ArrayList<Node> constructSolutionEAS(double[][] pm, double b, double a, double rho) {
		ArrayList<Integer> tour =  new ArrayList<>();
		int startIndex = generator.nextInt(nodes.size() + 1);
		int nextCity;

		tour.add(startIndex); //Initial node

		while(tour.size() != nodes.size()) {
			nextCity = pickNextCityEAS(tour, pm, b, a, rho, tour.get(tour.size()-1));
			tour.add(nextCity);
		}

		tour.add(startIndex);

		ArrayList<Node> tmp = new ArrayList<>();
		return tmp;
	}

	public static int pickNextCityEAS(ArrayList<Integer> tour , double[][] pmat , double b, double a, double rho, int lastCityID) {
		Node lastCity = nodes.get(tour.get(tour.size() - 1));
		ArrayList<Double> nodeProbs = new ArrayList<>();
		ArrayList<Integer> nodeTracker = new ArrayList<>();
		double denomSum = 0.0;

		//Calculate of all probs denominator and numerators of each node
		for(int i = 0; i < nodes.size(); i++){
			if(!tour.contains(i)){
				double invDistance = 1/(euclideanDistance2D(lastCity , nodes.get(i)));
				double pLevel = pmat[lastCityID][i];
				double nodeNum = (Math.pow(invDistance, b)) * (Math.pow(pLevel, a));
				// System.out.println("NodeNum: " + nodeNum);
				nodeProbs.add(nodeNum);
				nodeTracker.add(i);
				denomSum += nodeNum;
			}
		}

		for (int i = 0; i < nodeProbs.size(); i++){
			// System.out.println("Final: " + (nodeProbs.get(i)/denomSum));
			nodeProbs.set(i, (nodeProbs.get(i)/denomSum));
		}

		// System.exit(0);

		double cumulativeSum = 0;
		double random = generator.nextDouble();
		for (int i = 0; i < nodeProbs.size(); i++) {
			cumulativeSum += nodeProbs.get(i);
			if (random < cumulativeSum) {
				return nodeTracker.get(i);
			}
		}

		//Determine boundaries for each unvisited city

		// //Binary search to choose next city
		// double val = generator.nextDouble();
		// if(val <=nodeProbs.get(0)){
		// 	return 0;
		// } else{
		// 	int hi = nodeProbs.size();
		// 	int low = 0;
		// 	while(low <= hi){
		// 		int mid = low + (hi - low)/2;
		// 		if(val < nodeProbs.get(mid)){
		// 			hi = mid;
		// 		} else if(val > nodeProbs.get(mid)){
		// 			low = mid;
		// 		} else if((hi - low) == 1){
		// 			return nodeTracker.get(low);
		// 		}
		// 	}
		// }

		// System.out.println("Size: " + nodeProbs.size());

		// double sum = 0;
		// for (int i = 0; i < nodeProbs.size(); i++) {
		// 	sum+=nodeProbs.get(i);
		// }

		// System.out.println("Sum: " + sum);

		// System.out.println("HERE");

		return 0;

	}

	public static void globalPheromoneUpdateEAS(ArrayList<Node> candTour , ArrayList<Node> bestTour , double[][] phMatrix) {
		int costOne = computeCost(candTour);
		int costTwo = computeCost(bestTour);
		int k = 0;

		if(costOne == costTwo){ //MARCUS
			for (int i = 0; i < candTour.size(); i++){
				if(i ==(candTour.size() - 1)){
					k = 0;
				} else{
					k = k + 1;
				}
			}
		}

	}

	public static Boolean evaluateStopCondition(int stopCondition, int currIt, 
										long startTime, ArrayList<Node> tour) {

		switch(stopCondition) {
			case 0: 
				return (currIt < numberofIterations);
			case 1:	
				return (((computeCost(tour) / optimalTourCost) - 1) > errorAllowed);
			// case 2:
			// 	return ((currIt < numberofIterations) 
			// 		&& (((computeCost(tour) / optimalTourCost) - 1) > errorAllowed));
			case 2:
				return (((System.nanoTime() - startTime) / 1000000000.0) < secondsAllowed);
			case 3:
				return ((((computeCost(tour) / optimalTourCost) - 1) > errorAllowed)
					&& ((System.nanoTime() - startTime) / 1000000000.0) < secondsAllowed);
			default:
				printErrorAndExit();
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

	public static void initializeTauZero() {
		ArrayList<Integer> nearestNeighborTour = new ArrayList<>();
		int randomCity = generator.nextInt(nodes.size())+1;
		nearestNeighborTour.add(randomCity);

		while(nearestNeighborTour.size() != nodes.size()) {
			int closestNeighbor = getClosestCityTo(nearestNeighborTour.get(
										nearestNeighborTour.size()-1),nearestNeighborTour);
			nearestNeighborTour.add(closestNeighbor);
		}
		//convert Integer tour to Node tour.
		ArrayList<Node> rt = new ArrayList<>();
		for(int i = 0; i < nearestNeighborTour.size(); i++) {
			rt.add(getCity(nearestNeighborTour.get(i)));
		}
		tauZero = (1.0/(numberofAnts * computeCost(rt)));
	}

	public static int getClosestCityTo(int cityID, ArrayList<Integer> tour) {
		int closestCity = -1;
		double closestDistance = Integer.MAX_VALUE;
		for(int i = 0; i < nodes.size(); i++) {
			if(!tour.contains(i+1)) {
				double tempDist = euclideanDistance2D(getCity(cityID),getCity(i+1));
				if(tempDist < closestDistance) {
					closestDistance = tempDist;
					closestCity = (i+1);
				}
			}
		}
		// System.out.println("Closest City: " + closestCity);
		return closestCity;
	}

	public static Node getCity(int cityID) {
		for(int i = 0; i < nodes.size(); i++) {
			if(nodes.get(i).id == cityID) {
				return nodes.get(i);
			}
		}
		System.out.println("City: " + cityID + " not found. Exiting");
		System.exit(1);
		return null;
	}

	public static void processProblemLine(String[] tokens) {
		try {
			Integer.parseInt(tokens[0]); //make sure its a node line
			tokens = formatNodeInput(tokens);

			Node temp = new Node(Integer.parseInt(tokens[0]),
							Double.parseDouble(tokens[1]),
							Double.parseDouble(tokens[2]));
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