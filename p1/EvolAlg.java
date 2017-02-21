/*
MAXSAT - Project 1
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

*/

import java.io.*;
import java.util.*;

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

public class EvolAlg {

	//General Parameters
	private static String problemFilePath;
	private static Boolean whichAlgorithm; // 1 - GA, 0 - PBIL
	private static int numberOfVariables;
	private static int numberOfClauses;

	private static double mutationProbability;

	private static int bestIteration;
	private static int currIteration;

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

	private static long timeStart;
	private static long timeFinish;
	private static double timeElapsedSeconds;


	public static void main(String[] args) {
		if(args.length != 8) {
			System.out.println("Parameters were not supplied correctly");
			System.exit(1); //exit with error
		}
		readAndPrintParams(args);
		readFormula(problemFilePath);

		ArrayList<Double> sol = new ArrayList<>();

		if(whichAlgorithm) {
			//run GA
			ArrayList<Integer> temp;
			timeStart = System.nanoTime();
			temp = ga(numberOfIndividualsInThePopulation, breedingPoolSelectionMethod,
					crossoverMethod, crossoverProbability, mutationProbability,
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
						negativeLearningRate, mutationProbability, mutationAmount,
						numberOfIterations);
			timeFinish = System.nanoTime();
		}

		timeElapsedSeconds = (timeFinish - timeStart) / 1000000000.0;

		output(problemFilePath, numberOfVariables, numberOfClauses, -1, timeElapsedSeconds, sol);
	}

	public static int evaluateFitness(ArrayList<Integer> sample) {
		//Fitness = number of clauses met
		//Space = or
		//negative = not
		int fitness = 0;

		for (int i = 0; i < formula.size(); i++) { //For each clause
			for (int j = 0; j < formula.get(i).size(); j++) { //For each variable in clause
				if ((formula.get(i).get(j) > 0 && sample.get(Math.abs(formula.get(i).get(j))-1) == 1) 
						|| (formula.get(i).get(j) < 0 && sample.get(Math.abs(formula.get(i).get(j))-1) == 0)) {
					fitness++;
					break;
				}
			}
		}

		return fitness;
	}

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

	public static int bestSolutionIndex(ArrayList<Integer> fitnessEvaluations) {

		int bestSolutionIndex = 0;
		int currentMax = fitnessEvaluations.get(0);
		for(int i = 1; i < fitnessEvaluations.size(); i++) {
			if(fitnessEvaluations.get(i) > currentMax) {
				bestSolutionIndex = i;
				currentMax = fitnessEvaluations.get(i);
				bestIteration = currIteration;
			}
		}

		return bestSolutionIndex;

	}

	public static ArrayList<Integer> ga(int indInPopulation, String selection,
							String crossover, double crossoverProb, 
							double mutationProb, int numGenerations) {

		ArrayList<ArrayList<Integer>> population = initializePopulation(indInPopulation);

		ArrayList<Integer> fitnessEvaluations = new ArrayList<>();
		for(int i = 0; i < population.size(); i++) {
			fitnessEvaluations.add(evaluateFitness(population.get(i)));
		}

		ArrayList<Integer> best = population.get(bestSolutionIndex(fitnessEvaluations));

		ArrayList<ArrayList<Integer>> parents = new ArrayList<>();

		while(numGenerations > 0) {

			currIteration++;
			// System.out.println("Parents Size: " + parents.size());

			//Select parents
			if(selection.equals("rs")) { //rank selection
				parents = rs(population);
			} 
			else if(selection.equals("ts")) { //tournament selection
				parents = ts(population);
			}
			else if(selection.equals("bs")) { //boltzmann selection
				parents = bs(population);
			}
			else {
				System.exit(1);
			}

			// System.out.println("Parents Size: " + parents.size());

			// for (int i = 0; i < parents.size(); i++) {
			// 	System.out.println("Parent " + i + " Size: " + parents.get(i).size());
			// }
			// System.out.println();

			// System.out.println("Parents Size: " + parents.size());

			ArrayList<ArrayList<Integer>> children = new ArrayList<>();
			Child newChildren = new Child();

			for (int i = 0; i < parents.size(); i+=2) {
				if (crossover.equals("1c")) {
					newChildren = onepoint(parents.get(i), parents.get(i+1), crossoverProbability);
				} else if (crossover.equals("uc")) {
					newChildren = uniform(parents.get(i), parents.get(i+1), crossoverProbability);
				} else {
					System.exit(1);
				}

				ArrayList<Integer> c1 = newChildren.child1;
				ArrayList<Integer> c2 = newChildren.child2;

				c1 = mutateChild(c1);
				c2 = mutateChild(c2);

				children.add(c1); //Add children to children array list
				children.add(c2);

				// System.out.println(c1);
				// System.out.println(c2);

				// System.out.println("Children Size: " + children.size());

			}

			// System.out.println("Children Size: " + children.size());

			// for (int i = 0; i < children.size(); i++) {
			// 	System.out.println(children.get(i));
			// }


			ArrayList<Integer> childrenFitnessEvaluations = new ArrayList<>();
			for(int i = 0; i < children.size(); i++) {
				childrenFitnessEvaluations.add(evaluateFitness(children.get(i)));
				// System.out.println()
			}

			best = children.get(bestSolutionIndex(childrenFitnessEvaluations));

			// System.out.println("Fitness: " + evaluateFitness(best));

			population = children;

			// System.out.println("Population Size: " + population.size());

			numGenerations--;
		}

		return best;
	}

	public static ArrayList<Integer> mutateChild(ArrayList<Integer> child) {

		for (int i = 0; i < child.size(); i++) {
			if ((child.get(i) == 1) && (generator.nextDouble() < mutationProbability)) { //child at i is 1
				child.set(i, 0);
			}
			else if (generator.nextDouble() < mutationProbability) { //child at i is 0
				child.set(i, 1);
			}
		}

		// System.out.println("Child Size: " + child.size());

		return child;
	}

	public static Child onepoint(ArrayList<Integer> parent1, ArrayList<Integer> parent2, 
												double crossoverProbability) {
		Child c = new Child();

		if (generator.nextDouble() > crossoverProbability) {
			c.child1 = parent1;
			c.child2 = parent2;
			return c;
		}
		int randIndex = generator.nextInt(parent1.size());
		for(int i = 0; i < parent1.size(); i++){
			if(i >= randIndex){
				c.addToChild1(parent2.get(i));
				c.addToChild2(parent1.get(i));
			}
			else {
				c.addToChild1(parent1.get(i));
				c.addToChild2(parent2.get(i));
			}
		}

		// System.out.println("Child 1 Size: " + c.child1.size());
		// System.out.println("Child 2 Size: " + c.child2.size());

		return c;
	}

	//Check is cross over.
	//Equal prob pick from parent.
	public static Child uniform(ArrayList<Integer> parent1, ArrayList<Integer> parent2, 
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

	public static ArrayList<ArrayList<Integer>> rs(ArrayList<ArrayList<Integer>> population) {
		//Want selected.size == population.size
		ArrayList<ArrayList<Integer>> selected = new ArrayList<>();

		int totalRank = ((population.size()*(population.size() + 1))/2);

		// System.out.println("Total: " + totalRank);

		ArrayList<rankedIndividual> fitnessEvaluations = new ArrayList<>();
		for(int i = 0; i < population.size(); i++) {
			rankedIndividual newIndividual = new rankedIndividual();
			newIndividual.setIndividual(population.get(i));
			newIndividual.setFitness(evaluateFitness(population.get(i)));
			// System.out.println("New Fitness: " + newIndividual.getFitness());
			fitnessEvaluations.add(newIndividual);
			// System.out.println("Size: " + fitnessEvaluations.size());
			// System.out.println("Fitness of New: " + fitnessEvaluations.get(fitnessEvaluations.size()-1).getFitness());
		}

		// for (int i = 0; i < fitnessEvaluations.size(); i++) {
		// 	System.out.println("I: " + i);
	 //    	System.out.println("Fitness B4: " + fitnessEvaluations.get(i).getFitness());
	 //    }

		Collections.sort(fitnessEvaluations, new Comparator<rankedIndividual>() {
	        @Override public int compare(rankedIndividual i1, rankedIndividual i2) {
	            return i1.getFitness() - i2.getFitness();
	        }
	    });

	    // for (int i = 0; i < fitnessEvaluations.size(); i++) {
	    // 	System.out.println("Fitness: " + fitnessEvaluations.get(i).fitness);
	    // }

	    for (int i = 0; i < fitnessEvaluations.size(); i++) {
			fitnessEvaluations.get(i).setProbability(((double)(i+1) / totalRank));
	    }

	  //   for (int i = 0; i < fitnessEvaluations.size(); i++) {
			// System.out.println("Probability: " + fitnessEvaluations.get(i).getProbability());
	  //   }

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

	//M = 2, K = 1
	public static ArrayList<ArrayList<Integer>> ts(ArrayList<ArrayList<Integer>> population) {
		ArrayList<ArrayList<Integer>> selected = new ArrayList<>();

		while(selected.size() < population.size()) {
			ArrayList<Integer> option1 = population.get(generator.nextInt(numberOfVariables));
			ArrayList<Integer> option2 = population.get(generator.nextInt(numberOfVariables));

			int f1 = evaluateFitness(option1);
			int f2 = evaluateFitness(option2);

			selected.add(((f1 > f2) ? option1 : option2));
		}

		return selected;
	}

	public static ArrayList<ArrayList<Integer>> bs(ArrayList<ArrayList<Integer>> population) {
		ArrayList<ArrayList<Integer>> selected = new ArrayList<>();

		// System.out.println("Here1");


		// System.out.println("NUmber of Clauses: " + numberOfClauses);

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
			// System.out.println("Probability B4: " + prob);
			rankedIndividual newIndividual = new rankedIndividual();
			newIndividual.setIndividual(population.get(i));
			newIndividual.setProbability(prob);
			fitnessEvaluations.add(newIndividual);
			// System.out.println("Probability: " + newIndividual.getProbability());
		}

		// System.out.println("Fitness Ev");

		// for (int i = 0; i < fitnessEvaluations.size(); i++) {
		// 	System.out.println("I: " + i);
		// 	System.out.println("Probability: " + fitnessEvaluations.get(i).getProbability());
	 //    }

	    // System.out.println("Here6");

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

	    // System.out.println("Here4");

	  //   for (int i = 0; i < fitnessEvaluations.size(); i++) {
			// System.out.println("Probability Sorted: " + fitnessEvaluations.get(i).getProbability());
	  //   }

	    // System.out.println("Here5");

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

	public static ArrayList<Double> pbil(int indPerIteration, double posLearningRate, 
							double negLearningRate, double mutationProb,
							double mutationAmt, int numIterations) {

		ArrayList<Double> probVector = new ArrayList<>();
		for(int i = 0; i < numberOfVariables; i++) {
			probVector.add(0.5);
		}

		int mostFit = 0;

		int bestVectorIndex, worstVectorIndex;

		while(numIterations > 0) {

			currIteration++;

			//generate individuals
			ArrayList<ArrayList<Integer>> samples = new ArrayList<>();
			ArrayList<Integer> fitnessEvaluations = new ArrayList<>();

			for(int i = 0; i < indPerIteration; i++) { //for each individual
				ArrayList<Integer> sample = new ArrayList<>();
				for(int j = 0; j < numberOfVariables; j++) { //for each variable

					// double random = generator.nextDouble();
					// System.out.println("Random: " + random);

					sample.add(((generator.nextDouble() > probVector.get(j)) ? 0 : 1));
				}

				samples.add(sample);
				fitnessEvaluations.add(evaluateFitness(samples.get(i)));
				// System.out.println("Fitness: " + fitnessEvaluations.get(i));
			}

			// System.out.println("Samples Size: " + samples.size());

			bestVectorIndex = 0;
			worstVectorIndex = 0;

			for (int i = 1; i < fitnessEvaluations.size(); i++) {
				if (fitnessEvaluations.get(i) < fitnessEvaluations.get(worstVectorIndex)) {
					// System.out.println("Worse fitness: " + fitnessEvaluations.get(worstVectorIndex) + ">" + fitnessEvaluations.get(i));
					worstVectorIndex = i;
				}
				if (fitnessEvaluations.get(i) > fitnessEvaluations.get(bestVectorIndex)) {
					// System.out.println("Better fitness: " + fitnessEvaluations.get(bestVectorIndex) + "<" + fitnessEvaluations.get(i));
					bestVectorIndex = i;
				}
			}

			// System.out.println("Best Index: " + bestVectorIndex);
			// System.out.println("Worst Index: " + worstVectorIndex);

			// if (mostFit != fitnessEvaluations.get(bestVectorIndex)) {
			// 	mostFit = fitnessEvaluations.get(bestVectorIndex);
			// 	System.out.println("FITNESS CHANGE");
			// }

			// System.out.println("Most fit: " + fitnessEvaluations.get(bestVectorIndex));
			// System.out.println("Least fit: " + fitnessEvaluations.get(worstVectorIndex));

			// System.out.println("Prob Vector: ");
			// for(int i = 0; i < numberOfVariables; i++) {
			// 	System.out.println(probVector.get(i));
			// }

			for (int i = 0; i < numberOfVariables; i++) {
				probVector.set(i,(probVector.get(i) * (1.0 - posLearningRate) + 
								(samples.get(bestVectorIndex).get(i) * posLearningRate)));
			}

			// System.out.println("Prob Vector: ");
			// for(int i = 0; i < numberOfVariables; i++) {
			// 	System.out.println(probVector.get(i));
			// }

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

			for (int i = 0; i < numberOfVariables; i++) {
				if (probVector.get(i) > 0.9) {
						probVector.set(i,0.9);
				} 
				else if (probVector.get(i) < 0.1) {
					probVector.set(i,0.1);
				}	
			}

			// System.out.println("Prob Vector: ");
			// for(int i = 0; i < numberOfVariables; i++) {
			// 	System.out.println("Variable " + i + ": " + probVector.get(i));
			// }

			numIterations--;

			// ArrayList<Integer> processed = new ArrayList<>();

			// // System.out.println("Prob vector");

			// for(int i = 0; i < probVector.size(); i++) {

			// 	// System.out.println(probVector.get(i));

			// 	if(probVector.get(i) > 0.5) {
			// 		processed.add(1);
			// 	} else {
			// 		processed.add(0);
			// 	}
			// }

			// System.out.println();	
			// System.out.println();		

			// int satisfied = evaluateFitness(processed);
			// System.out.println("# of Satisfied Clauses: " + satisfied);

		}
		
		return probVector;

	}

	public static void output(String problemFP, int numVars, int numClauses, 
							  int iteration, double seconds, ArrayList<Double> sol) {
		// System.out.println("Job File: "+problemFP);
		// System.out.println("# of Variables: " + numVars);
		// System.out.println("# of Clauses: " + numClauses);

		//create the solution vector, i.e. process the probabilities, maybe
		//this should be done before calling output
		ArrayList<Integer> processed = new ArrayList<>(); //placeholder

		for(int i = 0; i < sol.size(); i++) {
			if(sol.get(i) > 0.5) {
				processed.add(1);
			} else {
				processed.add(0);
			}
		}

		int satisfied = evaluateFitness(processed);
		// System.out.println("# of Satisfied Clauses: " + satisfied);

		double percentage = (double)satisfied / (double)numClauses * 100;
		// System.out.println("% of Satisfied Clauses: " + percentage + "%\n");

		int lineCounter = 0;
		for(int i = 0; i < processed.size(); i++) {
			// System.out.print("v"+(i+1)+": "+processed.get(i)+"\t");
			if(lineCounter == 9) {
				//System.out.println(); //ten variables per line
				lineCounter = 0;
			} else {
				lineCounter++;
			}
		}
		// System.out.print("\n");
		//we currently don't support this statistic.
		// System.out.println("Iteration: " + iteration);

		System.out.print(numberOfClauses + " "+ satisfied + " " + percentage + " " + bestIteration + " " + seconds);
	}

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
