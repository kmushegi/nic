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
	public ArrayList<Double> location;
	public ArrayList<Double> velocity;
	public ArrayList<Integer> neighbors; //stored by particle number, i?

	public ArrayList<Double> personalBestLocation;
	public double personalBestValue;

	public ArrayList<Double> neighborhoodBestLocation;
	public double neighborhoodBestValue;

	public int dim;

	private static Random generator = new Random();

	Particle(int d) {
		this.location = new ArrayList<>();
		this.velocity = new ArrayList<>();
		this.personalBestLocation = new ArrayList<>();
		this.neighborhoodBestLocation = new ArrayList<>();
		this.dim = d;
	}

	public Particle initializeParticle(double minSpeed, double maxSpeed,
									double minLocation, double maxLocation,
									double bound) {
		for(int i = 0; i < this.dim; i++) {
			double x_i = bound * generator.nextDouble();

			while(x_i < minLocation || x_i > maxLocation) {
				x_i = bound * generator.nextDouble();
			}

			this.location.add(x_i);
			//this.velocity = ; //initialize velocity
			this.personalBestValue = -Double.MAX_VALUE;
			this.neighborhoodBestValue = -Double.MAX_VALUE;
		}
		personalBestLocation = location;
		return this;
	}

	public void setNeighborhood(ArrayList<Integer> n) {
		this.neighbors = n;
	}

	public int getDimension() {
		return dim;
	}
}