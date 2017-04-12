/*
ACO for TSP - Project 3
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Ant Colony Optimization for the Traveling Salesman Problem,
Project 3. This file contains functions used for outputting various information.
*/

import java.io.*;
import java.util.*;

public class Logger {

	public static void printStringArray(String[] s) {
		for(int i = 0; i < s.length; i++) {
			System.out.print(s[i] + " ");
		}
		System.out.print("\n");
	}

	public static void printNodeArrayList(ArrayList<Node> s) {
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

	public static void printErrorAndExit(String msg) {
		System.out.println("Error: " + msg);
		System.exit(1); //exit with error
	}
	
}