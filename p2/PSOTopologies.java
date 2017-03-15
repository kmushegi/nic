/*
Comparison of PSO Topologies - Project 2
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

The code in this file contains the implementation of A Comparison of 
Neighborhood Topologies in Particle Swarm Optimization, as part of Project 2.
*/

import java.io.*;
import java.util.*;
import java.util.List;
import java.util.ArrayList;

public class PSOTopologies {

	//PSO Constants
	private static final double CONSTRICTION_FACTOR = 0.7298;
	private static final double PHI1 = 2.05;
	private static final double PHI2 = 2.05;  
	private static final double BOUND = 400;

	//Topology - Neighborhood sizes
	private static final int ringNeighborSize = 2;
	private static final int randNeighborSize = 4;
	private static final int vnNeighborSize = 4;

	//Parameters in order of acceptance from CL
	//some parameters are fixed but can be specified on CL
	private static String whichTopology;
	private static int swarmSize;
	private static int numberOfIterations = 10000;
	private static String whichFunction;
	private static int functionDimensionality = 30;

	//PSO Bounds - set based on evaluation function
	private static double minVelocity;
	private static double maxVelocity;
	private static double minLocation;
	private static double maxLocation;

	private static ArrayList<Double> bestValues = new ArrayList<>();

	//Containers
	private static ArrayList<Particle> particles = new ArrayList<>();
	private static ArrayList<Integer> neighborhood;
	private static double bestGlobalLocation[];
	private static double bestGlobalValue = Double.MAX_VALUE;

	//Messages
	private static String incorrectParams = "Parameters were not supplied " 
											+ "correctly. Refer to man below:";
	private static String paramsOrder = "Param Order: gl/ri/vn/ra swarmSize "
											+ "numIterations rok/ack/ras " 
											+ "funcDimensionality";
	private static String swarmSizeError = "Choose swarm size of 16,30 or 49";

	//Util
	private static Random generator = new Random();

	//Timekeeping
	private static long timeStart;
	private static long timeFinish;
	private static double timeElapsedSeconds;

	public static void main(String[] args) {
		readParams(args);
		printParams();
		System.out.println();

		initializeBounds(whichFunction);
		initializeParticles();
		initializeTopology(whichTopology);

		timeStart = System.nanoTime();
		pso(whichFunction, whichTopology, particles, numberOfIterations);
		timeFinish = System.nanoTime();

		timeElapsedSeconds = (timeFinish - timeStart) / 1000000000.0; //ns to seconds
		System.out.println("Time Elapsed: " + timeElapsedSeconds + " seconds");
		System.out.println("Best Value: " + bestGlobalValue);
		System.out.println("Location: ");
		printDoubleArray(bestGlobalLocation);
		// System.out.print(timeElapsedSeconds + " " + bestGlobalValue + " ");
		// printBestValues();
		// System.out.println();

	}

	public static void pso(String function, String topology,
						ArrayList<Particle> particles, int iterations) {

		// bestValues.add(bestGlobalValue);

		while(iterations > 0) { //for each iteration
			iterations--;
			for(int i = 0; i < particles.size(); i++) { //for each particle
				for(int j = 0; j < particles.get(i).getDimension(); j++) { //for each dimension
					double pa = particles.get(i).personalBestLocation[j] //Calculate difference between current and personal best position
								- particles.get(i).location[j];

					double ga = particles.get(i).neighborhoodBestLocation[j] //Calculate difference between current and neighborhood best position
								- particles.get(i).location[j];

					double vi1 = (particles.get(i).velocity[j] //Calculate new velocity using velocity update function
							+ (generator.nextDouble() * PHI1 * pa)
							+ (generator.nextDouble() * PHI2 * ga));

					vi1 *= CONSTRICTION_FACTOR; //Constrict new velocity

					//update velocity & position
					particles.get(i).velocity[j] = vi1;
					particles.get(i).location[j] += particles.get(i).velocity[j];

				}

				//compute how good current position is
				double currPositionValue = eval(function,particles.get(i));

				//update personal and global bests if necessary
				if(currPositionValue <= particles.get(i).personalBestValue) {//Updating personal best
					particles.get(i).personalBestValue = currPositionValue;
					particles.get(i).personalBestLocation = Arrays.copyOf(
							particles.get(i).location, functionDimensionality);

					if(particles.get(i).personalBestValue <= bestGlobalValue) {//Updating global best
						bestGlobalValue = particles.get(i).personalBestValue;
						bestGlobalLocation = Arrays.copyOf(
							particles.get(i).location,functionDimensionality);
					}
				}
				
				//update neighborhood bests if necessary
				for(int nh = 0; nh < particles.get(i).neighbors.length; nh++) {
					Particle nbor = particles.get(particles.get(i).neighbors[nh]);
					double nborValue = eval(function, nbor);

					if(nborValue <= particles.get(i).neighborhoodBestValue) { //Updating neighboorhood best
						particles.get(i).neighborhoodBestLocation = nbor.location;
						particles.get(i).neighborhoodBestValue = nborValue;
					}
				}

				if(currPositionValue <= particles.get(i).neighborhoodBestValue) { //Setting particle's global best to self, if it is the best in its neighborhood
					particles.get(i).neighborhoodBestLocation = Arrays.copyOf(
						particles.get(i).location, functionDimensionality);
					particles.get(i).neighborhoodBestValue = currPositionValue;
				}
			}

			// if (iterations % 1000 == 0) {
			// 	bestValues.add(bestGlobalValue);
			// }
		}
	}

	//evaluate Particle p with function function.
	public static double eval(String function, Particle p) {
		double ev = 0.0;

		if(function.equals("rok")) {
			ev = evalRosenblock(p);
		} else if(function.equals("ack")) {
			ev = evalAckley(p);
		} else if(function.equals("ras")) {
			ev = evalRastrigin(p);
		} else {
			printErrorAndExit();
		}
		return ev;
	}
	
	//***************************************************************
	//each eval* function is based on the code provided by 			*
	//Professor Majercik in the PSO in-class lab and has been 		*
	//modified to support > 2 dimensions							*
	//***************************************************************

	//evaluate Particle p using Rosenblock evaluation. 
	public static double evalRosenblock(Particle p) {
		double ev = 0.0;
		double c_i, c_iplusone;

		for(int i = 0; i < p.getDimension() - 1; i++) {
			c_i = p.location[i];
			c_iplusone = p.location[(i+1)];
			ev += ((100.0 * (Math.pow((c_iplusone - Math.pow(c_i, 2.0)),2.0))) 
					+ Math.pow((c_i-1.0),2.0));
		}

		return ev;
	}

	//evaluate Particle p using Rastrigin evaluation.
	public static double evalRastrigin(Particle p) {
		double ev = 0.0;
		double c_i;

		for(int i = 0; i < p.getDimension(); i++) {
			c_i = p.location[i];
			ev += ((c_i * c_i) - 10.0*Math.cos(2.0*Math.PI*c_i) + 10.0);
		}

		return ev;
	}

	//evaluate Particle p using Ackley evaluation.
	public static double evalAckley(Particle p) {
		double ev = 0.0;

		double firstSum = 0.0;
		double secondSum = 0.0;

		double c_i;

		for(int i = 0; i < p.getDimension(); i++) {
			c_i = p.location[i];
			firstSum += Math.pow(c_i,2);
			secondSum += Math.cos(2.0*Math.PI*c_i);
		}

		ev =  -20.0 * Math.exp(-0.2 * Math.sqrt(firstSum/p.getDimension())) 
				- Math.exp(secondSum/p.getDimension()) + 20.0 + Math.E;
		return ev;
	}

	//initialize the specified topology
	public static void initializeTopology(String topology) {
		if(topology.equals("gl")) {
			initializeGlobalTopology();
		} else if(topology.equals("ri")) {
			initializeRingTopology();
		} else if(topology.equals("vn")) {
			initializeVonNeumannTopology();
		} else if(topology.equals("ra")) {
			initializeRandomTopology();
		} else {
			printErrorAndExit();
		}
	}

	//Every particle is neighbor with every other particle
	public static void initializeGlobalTopology() {
		int neighborhood[] = new int[swarmSize - 1];

		for(int i = 0; i < swarmSize; i++) {
			int counter = 0; //Tracks index of next neighbor in neighborhood
			for(int k = 0; k < swarmSize ; k++){
				if(i == k) { //cannot be your own neighbor
					continue;
				}
				neighborhood[counter] = k;
				counter++;
			}
			particles.get(i).setNeighborhood(neighborhood);		
		}
	}

	//Every particle is neighbors with next to it in particles array
	public static void initializeRingTopology() {
		int neighborhood[] = new int[ringNeighborSize];

		for(int i = 0; i < swarmSize; i++) { 
			if(i == 0) { //If first particle in array
				neighborhood[0] = (swarmSize - 1);
				neighborhood[1] = (i + 1);
			} else if (i == (swarmSize - 1)) { //If last particle in array
				neighborhood[0] = (i - 1);
				neighborhood[1] = 0;
			} else {
				neighborhood[0] = (i - 1);
				neighborhood[1] = (i + 1);
			}
			particles.get(i).setNeighborhood(neighborhood);		
		}
	}

	//Arranged in grid, particles neighbors are above,below, and to either side
	public static void initializeVonNeumannTopology() {
		int neighborhood[] = new int[vnNeighborSize];
		//Create nxn grid. Grid size is the ceiling of the squareroot of the size of the swarm.
		//This ensures that the square grid can fit all particles.
		int gridDimension = (int) Math.ceil(Math.sqrt(swarmSize));
		int[][] grid = new int[gridDimension][gridDimension];

		int particleIndex = 0; //populate the grid for von Neumann
		for (int i = 0; i < gridDimension; i++) {
			for (int j = 0; j < gridDimension; j++) {
				if (particleIndex > swarmSize - 1) { //If we have already put all particle indicies in the grid, set entry to -1 (i.e. empty)
					grid[i][j] = -1;
				} else { //Otherwise, set to particle index
					grid[i][j] = particleIndex;
					particleIndex++;
				}
			}
		}

		int colIndex, rowIndex;
		for (int i = 0; i < gridDimension; i++) {
			for (int j = 0; j < gridDimension; j++) {

				if (grid[i][j] == -1) { //If the current entry is -1 (i.e. empty), we don't need to set its neighbors.
					break;
				}

				//Finding the particle's left neighbor
				colIndex = j-1; //Set colIndex to particle entrie's immediate left
				//Check if the entry to the left is a valid entry
				while (colIndex < 0 || grid[i][colIndex] == -1) { //Continue to go to left if left entry is invalid (empty or out of bounds)
					if (colIndex < 0) { //If left is out of bounds of grid, warp around, and set colIndex to most right column
						colIndex = gridDimension - 1;
					} else { //If left entry is -1 (i.e. empty), continue to advance left.
						colIndex--;
					}
				}
				neighborhood[0] = grid[i][colIndex]; //Set particle's left neighbor to valid entry.

				//All below directions follow a similar logic to the left direction.

				//Finding the particle's right neighbor
				colIndex = j+1;
				while (colIndex > gridDimension - 1 || grid[i][colIndex] == -1) {
					if (colIndex > gridDimension - 1) {
						colIndex = 0;
					} else {
						colIndex++;
					}
				}
				neighborhood[1] = grid[i][colIndex];

				//Finding the particle's up neighbor
				rowIndex = i-1;
				while (rowIndex < 0 || grid[rowIndex][j] == -1) {
					if (rowIndex < 0) {
						rowIndex = gridDimension - 1;
					} else {
						rowIndex--;
					}
				}
				neighborhood[2] = grid[rowIndex][j];

				//Finding the particle's down neighbor
				rowIndex = i+1;
				while (rowIndex > gridDimension - 1 || grid[rowIndex][j] == -1) {
					if (rowIndex > gridDimension - 1) {
						rowIndex = 0;
					} else {
						rowIndex++;
					}
				}
				neighborhood[3] = grid[rowIndex][j];

				particles.get((i * gridDimension) + j).setNeighborhood(neighborhood);	
			}
		}
	}

	//Every particle is assigned four random neighbors
	public static void initializeRandomTopology() {
		int neighborhood[] = new int[randNeighborSize];

		for(int i = 0; i < particles.size(); i++) { //For each particles
			for(int neighborhoodIndex = 0; neighborhoodIndex < 
									randNeighborSize; neighborhoodIndex++) { //Assign new neighbors to each particle
				//set neighboor to a random index between 0 and # of particles-1 with replacement
				neighborhood[neighborhoodIndex] = generator.nextInt(particles.size()-1);
			}
			particles.get(i).setNeighborhood(neighborhood);
		}

	}	

	//initialize particles with given PSO parameters
	public static void initializeParticles() {
		for(int i = 0; i < swarmSize; i++) {
			Particle t = new Particle(functionDimensionality);
			t.initializeParticle(minVelocity, maxVelocity, 
								minLocation, maxLocation, BOUND);
			particles.add(t);
		}
	}

	//initialize bounds based on the eval function & project specs.
	public static void initializeBounds(String function) {
		if(function.equals("rok")) {
			minLocation = 15.0;
			maxLocation = 30.0;
			minVelocity = -2.0;
			maxVelocity = 2.0;
		} else if(function.equals("ack")) {
			minLocation = 16.0;
			maxLocation = 32.0;
			minVelocity = -2.0;
			maxVelocity = 4.0;
		} else if(function.equals("ras")) {
			minLocation = 2.56;
			maxLocation = 5.12;
			minVelocity = -2.0;
			maxVelocity = 4.0;
		} else {
			printErrorAndExit();
		}
	}

	public static void readParams(String[] args) {
		if(args.length != 3 && args.length != 5) {
			printErrorAndExit();
		} else if(args.length == 3) {
			whichTopology = args[0];
			swarmSize = Integer.parseInt(args[1]);
			whichFunction = args[2];
		} else if(args.length == 5) {
			whichTopology = args[0];
			swarmSize = Integer.parseInt(args[1]);
			numberOfIterations = Integer.parseInt(args[2]);
			whichFunction = args[3];
			functionDimensionality = Integer.parseInt(args[4]);
		}

		bestGlobalLocation = new double[functionDimensionality];
	}

	public static void printParams() {
		System.out.println("Topology: " + whichTopology
			+ "\nSwarm Size: " + swarmSize
			+ "\n# of Iterations: " + numberOfIterations
			+ "\nFunction: " + whichFunction
			+ "\nFunction Dimensionality: " + functionDimensionality);
	}

	//***********************
	//	Utilility functions	*
	//***********************

	public static int[] integerArrayListToIntArray(ArrayList<Integer> v) {
		int[] intArray = new int[v.size()];
		for (int i = 0; i < v.size(); i++) {
    		intArray[i] = v.get(i);
		}
		return intArray;
	}

	public static void printIntArray(int[] v) {
		int lineCounter = 0;
		for(int i = 0; i < v.length; i++) {
			System.out.print(v[i] + "\t");
			if(lineCounter == 4) {
				System.out.println(); //ten variables per line
				lineCounter = 0;
			} else {
				lineCounter++;
			}
		}
		System.out.print("\n");
	}

	public static void printDoubleArray(double[] v) {
		int lineCounter = 0;
		for(int i = 0; i < v.length; i++) {

			System.out.print(v[i] + "\t");
			if(lineCounter == 4) {
				System.out.println(); //ten variables per line
				lineCounter = 0;
			} else {
				lineCounter++;
			}
		}
		System.out.print("\n");
	}

	public static void printBestValues() {
		for (int i = 0; i < bestValues.size(); i++) {
			System.out.print(bestValues.get(i) + " ");
		}
	}

	public static void printErrorAndExit() {
		System.out.println(incorrectParams);
		System.out.println(paramsOrder);
		System.exit(1); //exit with error
	}

}