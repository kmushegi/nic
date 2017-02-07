/*
MAXSAT - Project 1
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

*/

import java.io.*;

public class EvolAlg {

	private static String problemFilePath;
	private static int numberOfIndividualsToGenerate;
	private static double positiveLearningRate;
	private static double negativeLearningRate;
	private static double mutationProbability;
	private static double mutationAmount;
	private static int numberOfIterations;
	private static Boolean whichAlgorithm; // 1 - GA, 0 - PBIL


	public static void main(String[] args) {
		if(args.length != 8) {
			System.out.println("Parameters were not supplied correctly");
			System.exit(1); //exit with error
		}
		readAndPrintParams(args);

	}

	public static void readAndPrintParams(String[] args) {
		problemFilePath = args[0];
		numberOfIndividualsToGenerate = Integer.parseInt(args[1]);
		positiveLearningRate = Double.parseDouble(args[2]);
		negativeLearningRate = Double.parseDouble(args[3]);
		mutationProbability = Double.parseDouble(args[4]);
		mutationAmount = Double.parseDouble(args[5]);
		numberOfIterations = Integer.parseInt(args[6]);

		System.out.println("Problem File Path: " + problemFilePath +
			"\n# of Individuals to Generate: " + numberOfIndividualsToGenerate +
			"\nPositive Learning Rate: " + positiveLearningRate +
			"\nNegative Learning Rate: " + negativeLearningRate +
			"\nMutation Probability: " + mutationProbability +
			"\nMutation Amount: " + mutationAmount +
			"\nNumber of Iterations: " + numberOfIterations);

		System.out.print("Algorithm: ");

		if(args[7].equals("g")) {
			System.out.println("GA");
			whichAlgorithm = true;
		} else if(args[7].equals("p")) {
			System.out.println("PBIL");
			whichAlgorithm = false;
		} else {
			System.out.println("Algorithm specified incorrectly");
			System.exit(1); //exit with error
		}
	}
}