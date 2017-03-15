/*
Comparison of PSO Topologies - Project 2
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

The code in this file contains the implementation of the Particle class.
*/


import java.io.*;
import java.util.*;

class Particle {
	public double location[];
	public double velocity[];
	public int neighbors[]; //Index of neighbors of particle

	public double personalBestLocation[];
	public double personalBestValue;

	public double neighborhoodBestLocation[];
	public double neighborhoodBestValue;

	public int dim; //Dimension of problem

	private static Random generator = new Random();

	Particle(int d) {
		this.location = new double[d];
		this.velocity = new double[d];
		this.personalBestLocation = new double[d];
		this.neighborhoodBestLocation = new double[d];
		this.dim = d;
	}

	public void initializeParticle(double minSpeed, double maxSpeed,
									double minLocation, double maxLocation,
									double bound) {
		for(int i = 0; i < this.dim; i++) {
			double x_i = bound * generator.nextDouble(); 

			while(x_i < minLocation || x_i > maxLocation) { //Generate a random initial location for each particle's dimension,
				x_i = bound * generator.nextDouble();		//ensuring that it is within a specific bound.
			}

			this.location[i] = x_i;
			this.velocity[i] = (minSpeed + generator.nextDouble() //Set random intial velocity, within bounds
									* (maxSpeed - minSpeed)); 
			this.personalBestValue = Double.MAX_VALUE;
			this.neighborhoodBestValue = Double.MAX_VALUE;
		}
		this.personalBestLocation = location;
	}

	public void setNeighborhood(int[] n) {
		this.neighbors = Arrays.copyOf(n,n.length);
	}

	public int getDimension() {
		return dim;
	}

}