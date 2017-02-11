/*
MAXSAT - Project 1
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

*/

import java.io.*;
import java.util.*;

public class EvolAlg {

	//General Parameters
	private static String problemFilePath;
	private static Boolean whichAlgorithm; // 1 - GA, 0 - PBIL

	private static double mutationProbability;

	//GA Parameters
	private static int numberOfIndividualsInThePopulation;
	private static String breedingPoolSelectionMethod;
	private static String crossoverMethod;
	private static double crossoverProbability;
	private static int numberOfGenerations;


	//PBIL Parameters
	private static int numberOfIndividualsToGenerate;
	private static double positiveLearningRate;
	private static double negativeLearningRate;
	private static double mutationAmount;
	private static int numberOfIterations;


	//Containers
	private static ArrayList<ArrayList<Integer>> formula = new ArrayList<>();


	public static void main(String[] args) {
		if(args.length != 8) {
			System.out.println("Parameters were not supplied correctly");
			System.exit(1); //exit with error
		}
		readAndPrintParams(args);
		readFormula(problemFilePath);
		printFormula();
	}

	public static void readAndPrintParams(String[] args) {
		problemFilePath = args[0];

		System.out.print("Algorithm: ");
		if(args[7].equals("g")) {
			System.out.println("GA");
			whichAlgorithm = true;

			numberOfIndividualsInThePopulation = Integer.parseInt(args[1]);
			breedingPoolSelectionMethod = args[2];
			crossoverMethod = args[3];
			crossoverProbability = Double.parseDouble(args[4]);
			mutationProbability = Double.parseDouble(args[5]);
			numberOfGenerations = Integer.parseInt(args[6]);

			System.out.println("Problem File Path: " + problemFilePath
				+ "\n# of Individuals In the Population: " + numberOfIndividualsInThePopulation
				+ "\nBreeding Pool Selection: " + breedingPoolSelectionMethod
				+ "\nCrossover Method: " + crossoverMethod
				+ "\nCrossover Probability: " + crossoverProbability
				+ "\nMutation Probability: " + mutationProbability
				+ "\nNumber of Generations: " + numberOfGenerations);

		} else if(args[7].equals("p")) {
			System.out.println("PBIL");
			whichAlgorithm = false;

			numberOfIndividualsToGenerate = Integer.parseInt(args[1]);
			positiveLearningRate = Double.parseDouble(args[2]);
			negativeLearningRate = Double.parseDouble(args[3]);
			mutationProbability = Double.parseDouble(args[4]);
			mutationAmount = Double.parseDouble(args[5]);
			numberOfIterations = Integer.parseInt(args[6]);

			System.out.println("Problem File Path: " + problemFilePath
				+ "\n# of Individuals to Generate: " + numberOfIndividualsToGenerate
				+ "\nPositive Learning Rate: " + positiveLearningRate
				+ "\nNegative Learning Rate: " + negativeLearningRate
				+ "\nMutation Probability: " + mutationProbability
				+ "\nMutation Amount: " + mutationAmount
				+ "\nNumber of Iterations: " + numberOfIterations);
			
		} else {
			System.out.println("Algorithm specified incorrectly");
			System.exit(1); //exit with error
		}
	}

	public static void readFormula(String problemFP) {
		try(BufferedReader br = new BufferedReader(new FileReader(problemFP))) {
			String line;
			while((line = br.readLine()) != null) {
				if(line.charAt(0) != 'c' && line.charAt(0) != 'p') {
					ArrayList<Integer> temp = new ArrayList<>();
					String[] tokens = line.split(" ");
					for(String token : tokens) {
						if(Integer.parseInt(token) == 0) {
							formula.add(temp);
							break;
						} else {
							temp.add(Integer.parseInt(token));
						}
					}
				}
			}
		} catch (IOException e) {
			System.err.format("IOException %s%n",e);
		}
	}

	public static void printFormula() {
		for(int i = 0; i < formula.size(); i++) {
			for(int j = 0; j < formula.get(i).size(); j++) {
				System.out.print(formula.get(i).get(j) + " ");
			}
			System.out.println();
		}
	}
}