/*
MAXSAT - Project 1
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

The code in this file contains the implementation of Genetic Algorithm and
Population Based Incremental Learning Algorithm, as part of Project 1.
*/

import java.io.*;
import java.util.*;

//*****************************************************
//					HELPER CLASSES					  *
//*****************************************************

class Child {
	public ArrayList<Integer> child1;
	public ArrayList<Integer> child2;

	Child() {
		child1 = new ArrayList<>();
		child2 = new ArrayList<>();
	}

	void addToChild1(int newInt) {
		child1.add(newInt);
	}

	void addToChild2(int newInt) {
		child2.add(newInt);
	}

}

class rankedIndividual {
	public ArrayList<Integer> individual;
	public int fitness;
	public double probability;

	void setIndividual(ArrayList<Integer> newIndividual) {
		individual = newIndividual;
	}

	void setFitness(int newFitness) {
		fitness = newFitness;
	}

	void setProbability(double newProbability) {
		probability = newProbability;
	}

	ArrayList<Integer> getIndividual() {
		return individual;
	}

	int getFitness() {
		return fitness;
	}

	double getProbability() {
		return probability;
	}

}

//*****************************************************
//					  MAIN CLASS					  *
//*****************************************************

public class EvolAlg {

	//General Parameters
	private static String problemFilePath; //location of .cnf file
	private static Boolean whichAlgorithm; // 1 - GA, 0 - PBIL
	private static int numberOfVariables;
	private static int numberOfClauses;
	private static double mutationProbability;

	private static int bestIteration = 0;
	private static int currIteration = 0;

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

	//Util
	private static Random generator = new Random();

	//Timekeeping
	private static long timeStart;
	private static long timeFinish;
	private static double timeElapsedSeconds;


	public static void main(String[] args) {
		if(args.length != 8) {
			System.out.println("Parameters were not supplied incorrectly");
			System.exit(1); //exit with error
		}
		readAndPrintParams(args);
		readFormula(problemFilePath);

		ArrayList<Double> sol = new ArrayList<>(); //solution list

		if(whichAlgorithm) { //run GA
			ArrayList<Integer> temp;

			timeStart = System.nanoTime();
			temp = ga(numberOfIndividualsInThePopulation, 
						breedingPoolSelectionMethod, crossoverMethod, 
						crossoverProbability, mutationProbability,
						numberOfGenerations);
			timeFinish = System.nanoTime();

			//maintain type consistency
			for(Integer i : temp) {
				sol.add(i.doubleValue());
			}
		} else {
			//run PBIL
			timeStart = System.nanoTime();
			sol = pbil(numberOfIndividualsToGenerate, positiveLearningRate,
						negativeLearningRate, mutationProbability, 
						mutationAmount, numberOfIterations);
			timeFinish = System.nanoTime();
		}

		timeElapsedSeconds = (timeFinish - timeStart) / 1000000000.0; //ns to seconds

		output(problemFilePath, numberOfVariables, numberOfClauses, 
							bestIteration, timeElapsedSeconds, sol);
	}

	//fitness = number of clauses satisfied
	public static int evaluateFitness(ArrayList<Integer> sample) {
		int fitness = 0;

		for (int i = 0; i < formula.size(); i++) { //For each clause
			for (int j = 0; j < formula.get(i).size(); j++) { //For each variable in clause
				if ((formula.get(i).get(j) > 0 //if clause is satisfied
						&& sample.get(Math.abs(formula.get(i).get(j))-1) == 1) 
						|| (formula.get(i).get(j) < 0 
						&& sample.get(Math.abs(formula.get(i).get(j))-1) == 0)) 
				{
					fitness++; //increase fitness value
					break;
				}
			}
		}

		return fitness;
	}

	//randomly generate population with 50% probability of being 1 or 0
	public static ArrayList<ArrayList<Integer>> initializePopulation(int inds) {
		ArrayList<ArrayList<Integer>> pp = new ArrayList<>();

		for(int i = 0; i < inds; i++) {
			ArrayList<Integer> temp = new ArrayList<>();
			for(int j = 0; j < numberOfVariables; j++) {
				temp.add(((generator.nextDouble() > 0.5) ? 0 : 1));
			}
			pp.add(temp);
		}

		return pp;
	}

	//compute the index of the best solution so far
	public static int bestSolutionIndex(ArrayList<Integer> fitnessEvaluations) {
		int bestSolutionIndex = 0;
		int currentMax = fitnessEvaluations.get(0);
		for(int i = 1; i < fitnessEvaluations.size(); i++) {
			if(fitnessEvaluations.get(i) > currentMax) {
				bestSolutionIndex = i;
				currentMax = fitnessEvaluations.get(i);
			}
		}

		return bestSolutionIndex;

	}

	//Genetic Algorithm Implementation, following pseudocode in the paper
	public static ArrayList<Integer> ga(int indInPopulation, String selection,
							String crossover, double crossoverProb, 
							double mutationProb, int numGenerations) {

		ArrayList<ArrayList<Integer>> population = initializePopulation(indInPopulation);

		ArrayList<Integer> fitnessEvaluations = new ArrayList<>();
		for(int i = 0; i < population.size(); i++) {
			fitnessEvaluations.add(evaluateFitness(population.get(i)));
		}

		ArrayList<Integer> best = population.get(bestSolutionIndex(fitnessEvaluations));
		int currentBestFitness = -1;

		ArrayList<ArrayList<Integer>> parents = new ArrayList<>();

		while(numGenerations > 0) {
			currIteration++;

			//Select parents
			if(selection.equals("rs")) { //rank selection
				parents = rs(population);
			} else if(selection.equals("ts")) { //tournament selection
				parents = ts(population);
			} else if(selection.equals("bs")) { //boltzmann selection
				parents = bs(population);
			} else {
				System.exit(1);
			}

			ArrayList<ArrayList<Integer>> children = new ArrayList<>();
			Child newChildren = new Child();

			for (int i = 0; i < parents.size(); i+=2) {
				if (crossover.equals("1c")) { //One point crossover
					newChildren = onepoint(parents.get(i), parents.get(i+1), 
												crossoverProbability);
				} else if (crossover.equals("uc")) { //Uniform crossover
					newChildren = uniform(parents.get(i), parents.get(i+1), 
												crossoverProbability);
				} else {
					System.exit(1);
				}

				ArrayList<Integer> c1 = newChildren.child1;
				ArrayList<Integer> c2 = newChildren.child2;

				c1 = mutateChild(c1); //Mutating children
				c2 = mutateChild(c2);

				children.add(c1); //Add children to children array list
				children.add(c2);

			}

			ArrayList<Integer> childrenFitnessEvaluations = new ArrayList<>();
			for(int i = 0; i < children.size(); i++) {
				childrenFitnessEvaluations.add(evaluateFitness(children.get(i)));
			}

			best = children.get(bestSolutionIndex(childrenFitnessEvaluations));
			int newBestFit = evaluateFitness(best);
			
			//Determining if the new best is better than the current best
			if(currentBestFitness < newBestFit){
				bestIteration = currIteration;
				currentBestFitness = newBestFit;
			}

			population = children;

			numGenerations--;
		}

		return best;
	}

	//mutate child with provided probability
	public static ArrayList<Integer> mutateChild(ArrayList<Integer> child) {

		for (int i = 0; i < child.size(); i++) {
			if ((child.get(i) == 1) && 
				(generator.nextDouble() < mutationProbability)) { //child at i is 1
				child.set(i, 0);
			}
			else if (generator.nextDouble() < mutationProbability) { //child at i is 0
				child.set(i, 1);
			}
		}

		return child;
	}

	//one-point crossover as specified in the project handout
	public static Child onepoint(ArrayList<Integer> parent1, 
								 ArrayList<Integer> parent2, 
								 double crossoverProbability) {
		Child c = new Child();

		if (generator.nextDouble() > crossoverProbability) {
			c.child1 = parent1;
			c.child2 = parent2;
			return c;
		}

		int randIndex = generator.nextInt(parent1.size());
		for(int i = 0; i < parent1.size(); i++){
			if(i >= randIndex){ //Crossing over from index randIndex onwards
				c.addToChild1(parent2.get(i));
				c.addToChild2(parent1.get(i));
			}
			else {
				c.addToChild1(parent1.get(i));
				c.addToChild2(parent2.get(i));
			}
		}

		return c;
	}

	//uniform crossover as specified in the project handout
	public static Child uniform(ArrayList<Integer> parent1, 
								ArrayList<Integer> parent2, 
								double crossoverProbability) {
		Child c = new Child();

		if (generator.nextDouble() > crossoverProbability) {
			c.child1 = parent1;
			c.child2 = parent2;
			return c;
		}

		for(int i = 0; i < parent1.size(); i++){
			if(generator.nextDouble() > 0.5){
				c.addToChild1(parent1.get(i));
			}
			else {
				c.addToChild1(parent2.get(i));
			}
		}

		for(int i = 0; i < parent1.size(); i++){
			if(generator.nextDouble() > 0.5){
				c.addToChild2(parent1.get(i));
			}
			else {
				c.addToChild2(parent2.get(i));
			}
		}

		return c;
	}

	//rank selection as specified in the project handout
	public static ArrayList<ArrayList<Integer>> rs(ArrayList<ArrayList<Integer>> population) {

		ArrayList<ArrayList<Integer>> selected = new ArrayList<>();

		int totalRank = ((population.size()*(population.size() + 1))/2);

		ArrayList<rankedIndividual> fitnessEvaluations = new ArrayList<>();
		for(int i = 0; i < population.size(); i++) {
			rankedIndividual newIndividual = new rankedIndividual();
			newIndividual.setIndividual(population.get(i));
			newIndividual.setFitness(evaluateFitness(population.get(i)));
			fitnessEvaluations.add(newIndividual);
		}

		//Soring the rankedIndividual array list based on their fitness
		Collections.sort(fitnessEvaluations, new Comparator<rankedIndividual>() {
	        @Override public int compare(rankedIndividual i1, rankedIndividual i2) {
	            return i1.getFitness() - i2.getFitness();
	        }
	    });

	    for (int i = 0; i < fitnessEvaluations.size(); i++) {
			fitnessEvaluations.get(i).setProbability(((double)(i+1) / totalRank));
	    }

	    while (selected.size() < population.size()) {
	    	double cumulativeSum = 0;
	    	double random = generator.nextDouble();
	    	for (int i = 0; i < fitnessEvaluations.size(); i++) {
	    		cumulativeSum += fitnessEvaluations.get(i).getProbability();
	    		if (random < cumulativeSum) {
	    			selected.add(fitnessEvaluations.get(i).getIndividual());
	    			break;
	    		}
	    	}
	    }

		return selected;

	}

	//tournament selection as specified in the project handout, M = 2, K = 1
	public static ArrayList<ArrayList<Integer>> ts(ArrayList<ArrayList<Integer>> population) {
		ArrayList<ArrayList<Integer>> selected = new ArrayList<>();

		while(selected.size() < population.size()) {
			ArrayList<Integer> option1 = population.get(generator.nextInt(population.size()));
			ArrayList<Integer> option2 = population.get(generator.nextInt(population.size()));

			int f1 = evaluateFitness(option1);
			int f2 = evaluateFitness(option2);

			selected.add(((f1 > f2) ? option1 : option2));
		}

		return selected;
	}

	//boltzmann selection as specified in the project handout
	public static ArrayList<ArrayList<Integer>> bs(ArrayList<ArrayList<Integer>> population) {
		ArrayList<ArrayList<Integer>> selected = new ArrayList<>();

		double denominator = 0;
		double adjustedFitness;
		for(int i = 0; i < population.size(); i++) {
			adjustedFitness = (double)evaluateFitness(population.get(i))/numberOfClauses;
			denominator += Math.exp(adjustedFitness);
		}

		ArrayList<rankedIndividual> fitnessEvaluations = new ArrayList<>();
		double numerator;
		double prob;
		for(int i = 0; i < population.size(); i++) {
			adjustedFitness = (double)evaluateFitness(population.get(i))/numberOfClauses;
			numerator = Math.pow(Math.E, adjustedFitness);
			prob = numerator/denominator;
			rankedIndividual newIndividual = new rankedIndividual();
			newIndividual.setIndividual(population.get(i));
			newIndividual.setProbability(prob);
			fitnessEvaluations.add(newIndividual);
		}

		//Soring the rankedIndividual array list based on their probability of selection
		Collections.sort(fitnessEvaluations, new Comparator<rankedIndividual>() {
	        @Override public int compare(rankedIndividual i1, rankedIndividual i2) {
	        	if (i1.getProbability()<i2.getProbability()) {
	        		return -1;
	        	}
	        	else if (i1.getProbability()>i2.getProbability()) {
	        		return 1;
	        	}
	        	return 0;
	        }
	    });

	    while (selected.size() < population.size()) {
	    	double cumulativeSum = 0;
	    	double random = generator.nextDouble();
	    	for (int i = 0; i < fitnessEvaluations.size(); i++) {
	    		cumulativeSum += fitnessEvaluations.get(i).getProbability();
	    		if (random < cumulativeSum) {
	    			selected.add(fitnessEvaluations.get(i).getIndividual());
	    			break;
	    		}
	    	}
	    }

		return selected;

	}

	//PBIL Algorithm Implementation, following pseudocode in the paper
	public static ArrayList<Double> pbil(int indPerIteration, double posLearningRate, 
							double negLearningRate, double mutationProb,
							double mutationAmt, int numIterations) {

		ArrayList<Double> probVector = new ArrayList<>();
		for(int i = 0; i < numberOfVariables; i++) {
			probVector.add(0.5);
		}

		int mostFit = 0;
		int bestVectorIndex, worstVectorIndex;
		int currentBestFitness = -1;

		while(numIterations > 0) {
			currIteration++;

			ArrayList<ArrayList<Integer>> samples = new ArrayList<>();
			ArrayList<Integer> fitnessEvaluations = new ArrayList<>();

			for(int i = 0; i < indPerIteration; i++) { //for each individual
				ArrayList<Integer> sample = new ArrayList<>();
				for(int j = 0; j < numberOfVariables; j++) { //for each variable
					//generate an individual w/ value T/F based on probability vector
					sample.add(((generator.nextDouble() > probVector.get(j)) ? 0 : 1));
				}

				samples.add(sample);
				fitnessEvaluations.add(evaluateFitness(samples.get(i)));
			}

			bestVectorIndex = 0;
			worstVectorIndex = 0;

			for (int i = 1; i < fitnessEvaluations.size(); i++) {
				if (fitnessEvaluations.get(i) < fitnessEvaluations.get(worstVectorIndex)) {
					worstVectorIndex = i;
				}
				if (fitnessEvaluations.get(i) > fitnessEvaluations.get(bestVectorIndex)) {
					bestVectorIndex = i;
				}
			}

			int newBestFit = fitnessEvaluations.get(bestVectorIndex);
			
			if(currentBestFitness < newBestFit) { //update best iteration as necessary
				bestIteration = currIteration;
				currentBestFitness = newBestFit;
			}

			//update the probability vector
			for (int i = 0; i < numberOfVariables; i++) {
				probVector.set(i,(probVector.get(i) * (1.0 - posLearningRate) + 
								(samples.get(bestVectorIndex).get(i) * posLearningRate)));
			}

			for (int i = 0; i < numberOfVariables; i++) {
				if (samples.get(bestVectorIndex).get(i) != samples.get(worstVectorIndex).get(i)) {
					probVector.set(i,(probVector.get(i) * (1.0 - negLearningRate) + 
								(samples.get(bestVectorIndex).get(i) * negLearningRate)));
				}
			}

			for (int i = 0; i < numberOfVariables; i++) {
				if (generator.nextDouble() < mutationProb) {
					int mutationDir = (generator.nextDouble() > 0.5) ? 1 : 0;
					probVector.set(i,(probVector.get(i) * (1.0 - mutationAmt) + 
								(mutationDir * mutationAmt)));
				}
			}


			//avoid extreme values in the probability vector
			for (int i = 0; i < numberOfVariables; i++) {
				if (probVector.get(i) > 0.9) {
						probVector.set(i,0.9);
				} 
				else if (probVector.get(i) < 0.1) {
					probVector.set(i,0.1);
				}	
			}

			numIterations--;

		}
		
		return probVector;

	}

	//output results in the format required by the project handout
	public static void output(String problemFP, int numVars, int numClauses, 
							  int iteration, double seconds, ArrayList<Double> sol) {

		System.out.println("Job File: "+problemFP);
		System.out.println("# of Variables: " + numVars);
 		System.out.println("# of Clauses: " + numClauses);

		//create the solution vector, i.e. process the probabilities
		ArrayList<Integer> processed = new ArrayList<>();

		for(int i = 0; i < sol.size(); i++) {
			if(sol.get(i) > 0.5) {
				processed.add(1);
			} else {
				processed.add(0);
			}
		}

		int satisfied = evaluateFitness(processed);
		System.out.println("# of Satisfied Clauses: " + satisfied);
		
		//compute percentage of satisfied clauses
		double percentage = (double)satisfied / (double)numClauses * 100;
		System.out.println("% of Satisfied Clauses: " + percentage + "%\n");

		int lineCounter = 0;
		for(int i = 0; i < processed.size(); i++) {
			System.out.print("v"+(i+1)+": "+processed.get(i)+"\t");
			if(lineCounter == 9) {
				System.out.println(); //ten variables per line
				lineCounter = 0;
			} else {
				lineCounter++;
			}
		}
		System.out.print("\n");
		System.out.println("Iteration: " + iteration);

		//testing print
		// System.out.print(numberOfClauses + " "+ satisfied + " " + percentage + " " + iteration + " " + seconds);
	}

	//parse command line arguments to the program based on specified algorithm
	public static void readAndPrintParams(String[] args) {
		problemFilePath = args[0];

		if(args[7].equals("g")) {
			whichAlgorithm = true;

			numberOfIndividualsInThePopulation = Integer.parseInt(args[1]);
			breedingPoolSelectionMethod = args[2];
			crossoverMethod = args[3];
			crossoverProbability = Double.parseDouble(args[4]);
			mutationProbability = Double.parseDouble(args[5]);
			numberOfGenerations = Integer.parseInt(args[6]);

			// System.out.println("Problem File Path: " + problemFilePath
			// 	+ "\n# of Individuals In the Population: " + numberOfIndividualsInThePopulation
			// 	+ "\nBreeding Pool Selection: " + breedingPoolSelectionMethod
			// 	+ "\nCrossover Method: " + crossoverMethod
			// 	+ "\nCrossover Probability: " + crossoverProbability
			// 	+ "\nMutation Probability: " + mutationProbability
			// 	+ "\nNumber of Generations: " + numberOfGenerations);

		} else if(args[7].equals("p")) {
			whichAlgorithm = false;

			numberOfIndividualsToGenerate = Integer.parseInt(args[1]);
			positiveLearningRate = Double.parseDouble(args[2]);
			negativeLearningRate = Double.parseDouble(args[3]);
			mutationProbability = Double.parseDouble(args[4]);
			mutationAmount = Double.parseDouble(args[5]);
			numberOfIterations = Integer.parseInt(args[6]);

			// System.out.println("Problem File Path: " + problemFilePath
			// 	+ "\n# of Individuals to Generate: " + numberOfIndividualsToGenerate
			// 	+ "\nPositive Learning Rate: " + positiveLearningRate
			// 	+ "\nNegative Learning Rate: " + negativeLearningRate
			// 	+ "\nMutation Probability: " + mutationProbability
			// 	+ "\nMutation Amount: " + mutationAmount
			// 	+ "\nNumber of Iterations: " + numberOfIterations);
			
		} else {
			System.out.println("Algorithm specified incorrectly");
			System.exit(1); //exit with error
		}
	}

	//Space = or | negative number = not (~)
	public static void readFormula(String problemFP) {
		try(BufferedReader br = new BufferedReader(new FileReader(problemFP))) {
			String line;
			while((line = br.readLine()) != null) {
				if(line.charAt(0) == 'p') {
					String[] tokens = line.split(" ");
					numberOfVariables = Integer.parseInt(tokens[2]);

					//some problems are formatted with an extra space
					if (tokens[3] != null && !tokens[3].isEmpty()) {
						numberOfClauses = Integer.parseInt(tokens[3]);
					}

					else {
						numberOfClauses = Integer.parseInt(tokens[4]);
					}

				} else if(line.charAt(0) != 'c' && line.charAt(0) != 'p') {
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
