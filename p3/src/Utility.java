/*
ACO for TSP - Project 3
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Ant Colony Optimization for the Traveling Salesman Problem,
Project 3. This file contains various utility functions used across the project.
*/

import java.io.*;
import java.util.*;

public class Utility {

	//some problem files contain trailing spaces; java tokenizer fails to recognize
	//this and spits out empty strings as tokens when splitting on " ". This method
	//detects such strings and does not include them in the output
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

	//compute euclidean distance between two nodes in the x,y plane
	public static double euclideanDistance2D(Node n1, Node n2) {
		return Math.sqrt(Math.pow(n1.x-n2.x,2) + Math.pow(n1.y-n2.y,2));
	}

	//compute cost of the given tour by adding up the distance between successive nodes
	public static int computeCost(ArrayList<Node> t) {
		int distance = 0;
		for(int i = 1; i < t.size(); i++) {
			Node n1 = t.get(i);
			Node n2 = t.get(i-1);
			distance += (euclideanDistance2D(n1,n2));
		}
		return distance;
	}

	//given a cityID return the corresponding Node object 
	public static Node getCity(int cityID, ArrayList<Node> nodes) {
		for(int i = 0; i < nodes.size(); i++) {
			if(nodes.get(i).id == cityID) {
				return nodes.get(i);
			}
		}
		System.out.println("City: " + cityID + " not found. Exiting");
		System.exit(1);
		return null;
	}

	//given a cityID return the ID of the closest city in the environment
	public static int getClosestCityTo(int cityID, ArrayList<Integer> tour, ArrayList<Node> nodes) {
		int closestCity = -1;
		double closestDistance = Integer.MAX_VALUE;
		for(int i = 0; i < nodes.size(); i++) {
			if(!tour.contains(i+1)) {
				double tempDist = euclideanDistance2D(getCity(cityID,nodes),getCity(i+1,nodes));
				if(tempDist < closestDistance) {
					closestDistance = tempDist;
					closestCity = (i+1);
				}
			}
		}
		return closestCity;
	}

	//convert a tour of cityIDs to a tour of Node objects
	public static ArrayList<Node> integerTourToNodeTour(ArrayList<Integer> t, ArrayList<Node> nodes) {
		ArrayList<Node> rt = new ArrayList<>();
		for(int i = 0; i < t.size(); i++) {
			rt.add(getCity(t.get(i),nodes));
		}
		return rt;
	}

	
}