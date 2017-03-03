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

public class PSOTopologies {

	private static final double CONSTRICTION_FACTOR = 0.7298;
	private static final double PHI1 = 2.05;
	private static final double PHI2 = 2.05;  
	private static final double BOUND = 400; //Generalized Schwefel 2.6
	
	private static final int ringNeighborSize = 2;
	private static final int randNeighborSize = 5;

	//Parameters in order of acceptance from CL
	//some parameters are fixed byt can be specified on CL
	private static String whichTopology;
	private static int swarmSize;
	private static int numberOfIterations = 10000;
	private static String whichFunction;
	private static int functionDimensionality = 30;

	//Params
	private static double minVelocity;
	private static double maxVelocity;
	private static double minLocation;
	private static double maxLocation;

	//Containers
	private static ArrayList<Particle> particles = new ArrayList<>();
	private static ArrayList<Integer> neighborhood;
	private static ArrayList<Double> bestL;
	private static double bestV;
	//Messages
	private static String incorrectParams = "Parameters were not supplied " 
											+ "correctly. Refer to man below:";
	private static String paramsOrder = "Param Order: gl/ri/vn/ra swarmSize "
											+ "numIterations rok/ack/ras " 
											+ "funcDimensionality";
	private static String swarmSizeError = "Choose swarm size of 16,30 or 49";

	//Util
	private static Random generator = new Random();

	public static void main(String[] args) {
		readParams(args);
		printParams();

		initializeBounds(whichFunction);
		initializeParticles();
		initializeTopology(whichTopology);

		PSO(whichFunction, whichTopology, particles, numberOfIterations);
	}

	public static void PSO(String function, String topology,
						ArrayList<Particle> particles, int iterations) {

		while(iterations > 0) {
			iterations--;

			for(int i = 0; i < particles.size(); i++) {

				Particle currentParticle = particles.get(i);
				for(int j = 0; j < currentParticle.getDimension(); j++) {
					double pa = currentParticle.personalBestLocation.get(j) 
								- currentParticle.location.get(j);

					double ga = currentParticle.neighborhoodBestLocation.get(j)
								- currentParticle.location.get(j);

					double vi1 = (currentParticle.velocity.get(j)
							+ (generator.nextDouble() * PHI1 * pa)
							+ (generator.nextDouble() * PHI2 * ga));
					vi1 *= CONSTRICTION_FACTOR;

					//update velocity
					particles.get(i).velocity.set(j,vi1);
					currentParticle.velocity.set(j,vi1);

					//update position based on velocity
					particles.get(i).location.set(j,(particles.get(i).location.get(j)+vi1));
					currentParticle.location.set(j,(particles.get(i).location.get(j)+vi1));
				}


				double currPositionValue = eval(function,currentParticle);

				if(currPositionValue < currentParticle.personalBestValue) {
					particles.get(i).personalBestValue = currPositionValue;
					particles.get(i).personalBestLocation = particles.get(i).location;
				}

				if(currentParticle.neighborhoodBestValue < bestV) {
					bestV = currentParticle.neighborhoodBestValue;
					bestL =currentParticle.neighborhoodBestLocation;
				}
			}

		}

	}

	public static  double eval(String function, Particle p) {
		double eval = 0.0;

		if(function.equals("rok")) {
			eval = evalRosenblock(p);
		} else if(function.equals("ack")) {
			eval = evalAckley(p);
		} else if(function.equals("ras")) {
			eval = evalRastrigin(p);
		} else {
			printErrorAndExit();
		}
		return eval;
	}

	public static double evalRosenblock(Particle p) {
		double ev = 0.0;

		for(int i = 0; i < p.getDimension() - 1; i++) {
			double c_i = p.location.get(i);
			double c_iplusone = p.location.get(i+1);
			ev += (100 * Math.pow((c_iplusone - c_i*c_i),2.0) 
					+ Math.pow((c_i-1.0),2.0));
		}

		return ev;
	}


	public static double evalRastrigin(Particle p) {
		double ev = 0.0;

		for(int i = 0; i < p.getDimension(); i++) {
			double c_i = p.location.get(i);
			ev += (c_i * c_i - 10.0*Math.cos(2.0*Math.PI*c_i) + 10.0);
		}

		return ev;
	}

	public static double evalAckley(Particle p) {
		double ev = 0.0;

		double firstSum = 0.0;
		double secondSum = 0.0;

		for(int i = 0; i < p.getDimension(); i++) {
			double c_i = p.location.get(i);

			firstSum += Math.pow(c_i,2);
			secondSum += Math.cos(2.0*Math.PI*c_i);
		}

		ev =  -20.0 * Math.exp(-0.2 * Math.sqrt(firstSum/2.0)) 
				- Math.exp(secondSum/2.0) + 20.0 + Math.E;
		return ev;
	}

	public static void initializeTopology(String topology) {
		if(topology.equals("gl")) {
			//init global topology, these should be separate functions
			initializeGlobalTopology();
		} else if(topology.equals("ri")) {
			//init ring topology
			initializeRingTopology();
		} else if(topology.equals("vn")) {
			//init von neumann topology
			initializeVonNeumannTopology();
		} else if(topology.equals("ra")) {
			//init random topology
			initializeRandomTopology();
		} else {
			printErrorAndExit();
		}
	}

	public static void initializeGlobalTopology() {
		neighborhood = new ArrayList<Integer>(swarmSize - 1);

		//For every particle
		for(int i = 0; i < swarmSize; i++){
			//There are swarmSize - 1 neighbors
			for(int k = 0; k < swarmSize - 1; k++){
				//Not neighbors with themselves
				if(i == k){
					continue;
				}
				neighborhood.set(k,i);
			}
			particles.get(i).setNeighborhood(neighborhood);
		}
	}

	public static void initializeRingTopology() {
		neighborhood = new ArrayList<Integer>(ringNeighborSize);

		//For every particle
		for(int i = 0; i < swarmSize; i++){
			if(i == 0){
				neighborhood.set(0 , (swarmSize - 1));
				neighborhood.set(1 , (i + 1));
			} else if (i == (swarmSize - 1)){
				neighborhood.set(0 , (i - 1));
				neighborhood.set(1 , 0);
			}else{
				neighborhood.set(0 , (i - 1));
				neighborhood.set(1 , (i + 1));
			}
			particles.get(i).setNeighborhood(neighborhood);		
			// System.out.print("Particle: "+ i + "\n"+"neighborhood: " + neighborhood);
		}
	}

	public static void initializeVonNeumannTopology() {
		
	}

	public static void initializeRandomTopology() {
		
	}

	public static void initializeParticles() {
		for(int i = 0; i < swarmSize; i++) {
			Particle t = new Particle(functionDimensionality);
			t.initializeParticle(minVelocity, maxVelocity, 
								minLocation, maxLocation, BOUND);
			particles.add(t);
		}
	}

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

		if(swarmSize != 16 && swarmSize != 30 && swarmSize != 49) {
			System.out.println(swarmSizeError);
			printErrorAndExit();
		}
	}

	public static void printErrorAndExit() {
		System.out.println(incorrectParams);
		System.out.println(paramsOrder);
		System.exit(1); //exit with error
	}

	public static void printParams() {
		System.out.println("Topology: " + whichTopology
			+ "\nSwarm Size: " + swarmSize
			+ "\n# of Iterations: " + numberOfIterations
			+ "\nFunction: " + whichFunction
			+ "\nFunction Dimensionality: " + functionDimensionality);
	}
}