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

class Particle {
	public int x;
	public int y;
	public int v;

	Particle() {
		x = -1;
		y = -1;
		v = 0;
	}

	Particle(int nx, int ny, int nv) {
		x = nx;
		y = ny;
		v = nv;
	}

	void setX(int nx) {
		x = nx;
	}

	void setY(int ny) {
		y = ny;
	}

	void setV(int nv) {
		v = nv;
	}

	int getX() {
		return x;
	}

	int getY() {
		return y;
	}

	int getV() {
		return v;
	}
}

public class PSOTopologies {

	//Parameters in order of acceptance from CL
	private static String whichTopology;
	private static int swarmSize;
	private static int numberOfIterations = 10000; //fixed but can be specifed on CL
	private static String whichFunction;
	private static int functionDimensionality = 30; //fixed but can be specifed on CL

	//Messages
	private static String incorrectParams = "Parameters were not supplied " 
											+ "correctly. Refer to man below:";
	private static String paramsOrder = "Param Order: gl/ri/vn/ra swarmSize "
											+ "numIterations rok/ack/ras " 
											+ "funcDimensionality";

	public static void main(String[] args) {
		readParams(args);
		printParams();
	}

	public static void readParams(String[] args) {
		if(args.length != 3 && args.length != 5) {
			System.out.println(incorrectParams);
			System.out.println(paramsOrder);
			System.exit(1); //exit with error
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
	}

	public static void printParams() {
		System.out.println("Topology: " + whichTopology
			+ "\nSwarm Size: " + swarmSize
			+ "\n# of Iterations: " + numberOfIterations
			+ "\nFunction: " + whichFunction
			+ "\nFunction Dimensionality: " + functionDimensionality);
	}
}