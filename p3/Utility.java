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

public class Utility {

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

	public static double euclideanDistance2D(Node n1, Node n2) {
		return Math.sqrt(Math.pow(n1.x-n2.x,2) + Math.pow(n1.y-n2.y,2));
	}

	public static int computeCost(ArrayList<Node> t) {
		int distance = 0;
		for(int i = 1; i < t.size(); i++) {
			Node n1 = t.get(i);
			Node n2 = t.get(i-1);
			distance += (euclideanDistance2D(n1,n2));
		}
		return distance;
	}

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

	public static ArrayList<Node> integerTourToNodeTour(ArrayList<Integer> t, ArrayList<Node> nodes) {
		ArrayList<Node> rt = new ArrayList<>();
		for(int i = 0; i < t.size(); i++) {
			rt.add(getCity(t.get(i),nodes));
		}
		return rt;
	}

	
}