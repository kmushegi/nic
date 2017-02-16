/*
MAXSAT - Project 1
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

*/
//comment


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
	private static int numberOfVariables;
	private static int numberOfClauses;


	//PBIL Parameters
	private static int numberOfIndividualsToGenerate;
	private static double positiveLearningRate;
	private static double negativeLearningRate;
	private static double mutationAmount;
	private static int numberOfIterations;


	//Containers
	private static ArrayList<ArrayList<Integer>> formula = new ArrayList<>();
	private static String bitString = "";


	public static void main(String[] args) {
		if(args.length != 8) {
			System.out.println("Parameters were not supplied correctly");
			System.exit(1); //exit with error
		}
		readAndPrintParams(args);
		readFormula(problemFilePath);
		printFormula();
		System.out.println(bitString);

		//IDEA: supply parameters into each algorithm so that global variables are only
		//used in the function call in order to avoid confusion
		//TODO: get feedback on this from Marcus & Ernesto
		if(whichAlgorithm) {
			//run GA
		} else {
			//run PBIL
		}
	}

	public static int evaluateFitness(ArrayList<Double> sample) {
		//Fitness = number of clauses met
		//Space = or
		//negative = not

		int fitness = 0;

		for (int i = 0; i < formula.size(); i++) { //For each clause
			for (int j = 0; j < formula[i].length; j++) { //For each variable in clause
				if ((formula[i][j] > 0 && sample[formula[i][j]] == 1) || (formula[i][j] < 0 && sample[formula[i][j]] == 0)) {
					fitness++;
					break;
				}
			}
		}

		return fitness;
	}

	public static ArrayList<Double> pbil(int indPerIteration, double posLearningRate, 
							double negLearningRate, double mutationProb,
							double mutationAmt, int numIterations) {

		ArrayList<Double> probVector = new ArrayList<>();
		for(int i = 0; i < numberOfVariables; i++) {
			probVector.add(0.5);
		}
		
		ArrayList<ArrayList<Integer>> samples = new ArrayList<>(indPerIteration);
		ArrayList<Double> fitnessEvaluations = new ArrayList<>(); //The higher the value, the better the fitness.

		while(numIterations > 0) {
			//generate individuals
			for(int i = 0; i < indPerIteration; i++) { //for each individual
				for(int j = 0; j < numberOfVariables; j++) { //for each variable
					samples[i][j] = (Math.random() > probVector[j]) ? 0 : 1;
				}
				fitnessEvaluations[i] = evaluateFitness(samples[i]);
			}

			int bestVectorIndex, worstVectorIndex;
			bestVectorIndex = worstVectorIndex = 0;

			for (int i = 1; i < fitnessEvaluations.size(); i++) {
				if (fitnessEvaluations[i] < fitnessEvaluations[worstVectorIndex]) {
					worstVectorIndex = i;
				}
				if (fitnessEvaluations[i] > fitnessEvaluations[bestVectorIndex]) {
					bestVectorIndex = i;
				}
			}

			for (int i = 0; i < numberOfVariables; i++) {
				probVector[i] = probVector[i] * (1.0 - posLearningRate) + (samples[bestVectorIndex][i] * posLearningRate);
			}

			for (int i = 0; i < numberOfVariables; i++) {
				if (samples[bestVectorIndex][i] != samples[worstVectorIndex][i]) {
					probVector[i] = probVector[i] * (1.0 - negLearningRate) + (samples[bestVectorIndex][i] * negLearningRate);
				}
			}

			for (int i = 0; i < numberOfVariables; i++) {
				if (Math.random() < mutationProb) {
					int mutationDir = (Math.random() > 0.5) ? 1 : 0;
					probVector[i] = probVector[i] * (1.0 - mutationAmt) + (mutationDir * mutationAmt);
				}
			}

			numIterations--;

		}
		
		return probVector;

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
				if(line.charAt(0) == 'p') {
					String[] tokens = line.split(" ");
					numberOfVariables = Integer.parseInt(tokens[2]);
					numberOfClauses = Integer.parseInt(tokens[3]);
				} else if(line.charAt(0) != 'c' && line.charAt(0) != 'p') {
					ArrayList<Integer> temp = new ArrayList<>();
					String[] tokens = line.split(" ");
					for(String token : tokens) {
						if(Integer.parseInt(token) == 0) {
							formula.add(temp);
							break;
						} else {
							// if(Integer.parseInt(token) < 0) {
							// 	bitString += "0";
							// } else {
							// 	bitString += "1";
							// }
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
