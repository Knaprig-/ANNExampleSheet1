import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

@SuppressWarnings("ALL")
public class MainNetwork {
	
	private Random rand = new Random();
	
	//Exercise 1
	private final int N = 200;
	private final int P = 40;
	private double[][] w = new double[N][N];
	private int[][] patterns = new int[P][N];
	private int[] network = new int[N];
	
	//Exercise 2
	private final double BETA2 = 2;
	
	//Exercise 3
	private final double ETA = 0.02;
	private final double BETA3 = 0.5;
	private double[] inputs3a = new double[2];
	private double output3a;
	private double[] weights3a = new double[2];
	private double threshold3a = 0;
	private LinkedList<Truple> trainingSet;
	private LinkedList<Truple> validationSet;
	private LinkedList<Tuple> trainingEnergyList = new LinkedList<>();
	private LinkedList<Tuple> validationEnergyList = new LinkedList<>();

	
	/**
	 * A simple two-value tuple.
	 */
	private class Tuple {
		public double x;
		public double y;
		
		public Tuple(double x, double y) {
			this.x = x;
			this.y = y;
		}
	}
	
	/**
	 * A three-value tuple, or a truple.
	 */
	private class Truple {
		public double x;
		public double y;
		public double z;
		
		public Truple(double x, double y, double z) {
			this.x = x;
			this.y = y;
			this.z = z;
		}
	}
	
	/**
	 * Signum function. Returns 1 when argument is positive or 0, -1 otherwise.
	 */
	public int signum(double x) {
		if (x < 0) {
			return -1;
		} else {
			return +1;
		}
	}
	
	//-------------Exercise 1----------------
	
	/**
	 * EXERCISE 1
	 * Performs task 1b.
	 */
	public void task1b() {
		int faultyUpdates = 0;
		for (int i = 0; i < 500; i++) {
			//Randomly generated 50/50
			generatePatterns();
			
			//Hebb's Rule
			calculateWeights();
			
			//Feed a given pattern.
			feed(1);
			
			synchronousDeterministicUpdate(1);
			
			for (int j = 0; j < N; j++) {
				if (patterns[0][j] != network[j]) {
					faultyUpdates++;
				}
			}
			if (i % 25 == 0) {
				System.out.print(".");
			}
		}
		System.out.println((double) faultyUpdates / 100000);
	}
	
	/**
	 * EXERCISE 1
	 * Generates random patterns,
	 * consisting of bits that are either -1 or +1.
	 */
	private void generatePatterns() {
		for (int i = 0; i < patterns.length; i++) {
			for (int j = 0; j < patterns[0].length; j++) {
				patterns[i][j] = rand.nextBoolean() ? 1 : -1;
			}
		}
	}
	
	/**
	 * EXERCISE 1 & 2
	 * Uses Hebb's rule to calculate the values of the weight matrix.
	 */
	private void calculateWeights() {
		for (int i = 0; i < w.length; i++) {
			for (int j = 0; j < w.length; j++) {
				double sum = 0;
				for (int k = 0; k < P; k++) {
					sum += patterns[k][i] * patterns[k][j];
				}
				sum /= N;
				w[i][j] = sum;
				//Comment this out if running exercise 1:
				if (i == j) {
					w[i][j] = 0;
				}
			}
		}
	}
	
	/**
	 * EXERCISE 1
	 * Sets the network's neurons to the same
	 * value as those in a stored pattern.
	 *
	 * @param index Index of the pattern to be fed
	 */
	private void feed(int index) {
		//arrays indexed by 0;
		network = new int[N];
		index--;
		for (int i = 0; i < N; i++) {
			network[i] = patterns[index][i];
		}
	}
	
	/**
	 * EXERCISE 1
	 * Updates the network synchronously for the given amount of steps.
	 */
	private void synchronousDeterministicUpdate(int steps) {
		for (int step = 0; step < steps; step++) {
			int[] temp = new int[N];
			
			for (int i = 0; i < N; i++) {
				double sum = 0;
				for (int j = 0; j < N; j++) {
					sum += w[i][j] * network[j];
				}
				int result = signum(sum);
				temp[i] = result;
			}
			
			network = temp;
		}
	}
	
	/**
	 * DEBUG
	 */
	private void printPatterns() {
		for (int i = 0; i < P; i++) {
			for (int j = 0; j < N; j++) {
				if (patterns[i][j] == 1) {
					System.out.print("[+]");
				} else {
					System.out.print("[-]");
				}
			}
			System.out.println("");
		}
	}
	
	/**
	 * DEBUG
	 */
	private void printNetwork() {
		System.out.println("Network state:");
		for (int i = 0; i < N; i++) {
			if (network[i] == 1) {
				System.out.print("[+]");
			} else {
				System.out.print("[-]");
			}
		}
		System.out.println("");
	}
	
	//---------------------Exercise 2-------------------------
	
	/**
	 * EXERCISE 2
	 * Performs tasks 2a and 2b.
	 */
	private void task2() throws IOException {
		final int ITERATIONS = 20;
		LinkedList<Tuple> steadyList = new LinkedList<>();
		boolean steadyFound = false;
		double maxX = 0;
		StringBuilder sb;
		
		Files.write(Paths.get("2aOutput.txt"), "".getBytes());
		
		for (int i = 1; i <= ITERATIONS; i++) {
			LinkedList<Tuple> mList = new LinkedList<>();
			sb = new StringBuilder();
			generatePatterns();
			calculateWeights();
			feed(1);
			int step = -1;
			double m = 0;
			double averageM = 0;
			int BATCH_SIZE = 50;
			steadyFound = false;
			
			while (true) {
				step++;
				asynchronousStochasticUpdate(1);
				m += getM(0);
				
				//When we have calculated a certain number of m values, it's time to average them
				if (step % BATCH_SIZE == 0) {
					//Don't divide m if we're on the first step
					if (step != 0) {
						m /= BATCH_SIZE;
					}
					averageM = ((averageM * mList.size()) + (m)) / (mList.size() + 1);
					mList.add(new Tuple(step, averageM));
					
					//Compares the latest value of the order parameter with the average of the six before it
					if (mList.size() > 7 && !steadyFound) {
						if (Math.abs(
								mList.get(mList.size() - 1).y -
										(mList.get(mList.size() - 2).y +
												mList.get(mList.size() - 3).y +
												mList.get(mList.size() - 4).y +
												mList.get(mList.size() - 5).y +
												mList.get(mList.size() - 6).y +
												mList.get(mList.size() - 7).y) / 6) < 0.0001) {
							
							steadyList.add(new Tuple(step, 0));
							steadyFound = true;
						}
					}
					//Forcibly exits if we exceed 10^5 iterations
					if (step == 100000) {
						maxX = 100000;
						break;
					}
					m = 0;
				}
			}
			
			//Build the string used to plot our data in Mathematica
			sb.append("g" + i + " = ListPlot[{");
			for (int j = 0; j < mList.size(); j++) {
				sb.append("{" + mList.get(j).x + "," + mList.get(j).y + "}");
				if (j != mList.size() - 1) {
					sb.append(", ");
				}
			}
			sb.append("}, Joined -> True, PlotRange->{0,1}, PlotStyle->RGBColor[" +
					rand.nextDouble() + "," +
					rand.nextDouble() + "," +
					rand.nextDouble() + "]];\n\n");
			Files.write(Paths.get("2aOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);
		}
		
		double steadySteps = 0;
		for (Tuple t : steadyList) {
			System.out.println(t.x);
			steadySteps += t.x;
		}
		steadySteps /= steadyList.size();
		System.out.println("Average steady state steps: " + steadySteps);
		
		//More output for Mathematica
		sb = new StringBuilder();
		sb.append("Show[");
		for (int i = 1; i <= ITERATIONS; i++) {
			sb.append("g" + i + ",");
		}
		sb.append(" PlotRange->{{0, " + maxX + "},{0,1}}]");
		Files.write(Paths.get("2aOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);
	}
	
	/**
	 * EXERCISE 2
	 * Updates a single neuron, chosen at random.
	 *
	 * @param steps The number of iterations to be performed
	 */
	private void asynchronousStochasticUpdate(int steps) {
		for (int step = 0; step < steps; step++) {
			//Choose a random bit to update
			int i = rand.nextInt(N);
			if (rand.nextDouble() < activation2(i)) {
				network[i] = 1;
			} else {
				network[i] = -1;
			}
		}
	}
	
	/**
	 * EXERCISE 2
	 * Activation function.
	 */
	private double activation2(int i) {
		return 1 / (1 + Math.exp(-2 * BETA2 * b2(i)));
	}
	
	/**
	 * EXERCISE 2
	 * Calculates the local field for neuron i.
	 */
	private double b2(int i) {
		double sum = 0;
		for (int j = 0; j < N; j++) {
			sum += w[i][j] * network[j];
		}
		return sum;
	}
	
	/**
	 * EXERCISE 2
	 * Calculates the order parameter with respect to the given pattern.
	 */
	private double getM(int pattern) {
		double result = 0;
		for (int i = 0; i < N; i++) {
			result += network[i] * patterns[pattern][i];
		}
		result /= N;
		return result;
	}
	
	//--------------------Exercise 3--------------------
	
	/**
	 * EXERCISE 3
	 * Performs task 3a.
	 */
	private void task3a() throws IOException {
		//Prepare our data
		trainingSet = fileToList("training.txt");
		validationSet = fileToList("validation.txt");
		normalizeInputs(trainingSet);
		normalizeInputs(validationSet);
		double[] classificationsTraining = new double[10];
		double[] classificationsValidation = new double[10];

		//Reset exported data.
		Files.write(Paths.get("3aOutput.txt"), "".getBytes());

		for (int run = 1; run <= 10; run++) {
			trainingEnergyList.clear();
			validationEnergyList.clear();
			System.out.println(".");
			//Set up weights and threshold randomly
			for (int i = 0; i < 2; i++) {
				weights3a[i] = rand.nextDouble() * 0.4 - 0.2;
			}
			threshold3a = rand.nextDouble() * 2 - 1;
			
			//Run the network for a 10^6 iterations
			for (int i = 0; i < 1000 * 1000; i++) {
				int randNum = rand.nextInt(trainingSet.size());
				feedPattern3a(randNum, trainingSet);
				output3a = activation3a(b3a());
				double error = (trainingSet.get(randNum).z-output3a)*activation3aPrime(b3a());
				//Update weights
				weights3a[0] += ETA*error*inputs3a[0];
				weights3a[1] += ETA*error*inputs3a[1];
				//Update Threshold
				threshold3a += -ETA*error;
				if (i % 5000 == 0){
					double trainingEnergy = getEnergy(trainingSet);
					trainingEnergyList.add(new Tuple(i, trainingEnergy));
					double validationEnergy = getEnergy(validationSet);
					validationEnergyList.add(new Tuple(i, validationEnergy));
				}
			}
			exportEnergyData(trainingEnergyList, rand.nextDouble() + "," + rand.nextDouble() + ", 0.5", (2*run-1),"3aOutput.txt");
			exportEnergyData(validationEnergyList, "0.5," + rand.nextDouble() + "," + rand.nextDouble(), (2*run), "3aOutput.txt");
			classificationsTraining[run-1] = getC(trainingSet);
			classificationsValidation[run-1] = getC(validationSet);
		}
		//More output for Mathematica
		StringBuilder sb = new StringBuilder();
		sb.append("Show[");
		for (int i = 1; i <= 10; i++) {
			sb.append("g" + (2*i-1) + ",");
		}
		sb.append(" PlotRange->{{0, 1000000},{50,200}}]\n");
		Files.write(Paths.get("3aOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);
		sb = new StringBuilder();
		sb.append("Show[");
		for (int i = 1; i <= 10; i++) {
			sb.append("g" + (2*i) + ",");
		}
		sb.append(" PlotRange->{{0, 1000000},{50,200}}]");
		Files.write(Paths.get("3aOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);

		double averageClass = 0;
		double minimumClass = 400000;
		double varianceClass = 0;
		for (int i = 0; i < 10; i++) {
			averageClass += classificationsTraining[i];
			if (classificationsTraining[i] < minimumClass){
				minimumClass = classificationsTraining[i];
			}
			varianceClass += Math.pow(classificationsTraining[i],2);
		}
		averageClass = (averageClass / 10);
		varianceClass = (varianceClass/10) - Math.pow(averageClass,2);
		System.out.println("Trainingset values:\n Avr:" + averageClass + "\nmin: " + minimumClass + "\nvar:" + varianceClass);

		averageClass = 0;
		minimumClass = 400000;
		varianceClass = 0;
		for (int i = 0; i < 10; i++) {
			averageClass += classificationsValidation[i];
			if (classificationsValidation[i] < minimumClass){
				minimumClass = classificationsValidation[i];
			}
			varianceClass += Math.pow(classificationsValidation[i],2);
		}
		averageClass = (averageClass / 10);
		varianceClass = (varianceClass/10) - Math.pow(averageClass,2);
		System.out.println("Validationset values:\n Avr:" + averageClass + "\nmin: " + minimumClass + "\nvar:" + varianceClass);
	}

	private double getEnergy(LinkedList<Truple> list) {
		double result = 0;
		for (int i = 0; i < list.size(); i++) {
			feedPattern3a(i, list);
			output3a = activation3a(b3a());
			result += Math.pow(list.get(i).z - output3a, 2);
		}
		result /= 2;
		return result;
	}

	private void exportEnergyData(LinkedList<Tuple> list, String s, int index, String filePath) throws IOException {

		StringBuilder sb = new StringBuilder();

		sb.append("g" + index + " = ListPlot[{");
		for (int j = 0; j < list.size(); j++) {
			sb.append("{" + list.get(j).x + "," + list.get(j).y + "}");
			if (j != list.size() - 1) {
				sb.append(", ");
			}
		}
		sb.append("}, Joined -> True, PlotRange->{0,200}, PlotStyle->RGBColor[" +
				s + "]];\n\n");
		Files.write(Paths.get(filePath), sb.toString().getBytes(), StandardOpenOption.APPEND);
	}

	private double getC(LinkedList<Truple> dataSet) {
		double result = 0;
		for (int i = 0; i < dataSet.size(); i++) {
			feedPattern3a(i, dataSet);
			output3a = activation3a(b3a());
			result += Math.abs(signum(output3a) - dataSet.get(i).z);
		}
		result /= 2*dataSet.size();
		return result;
	}

	private double getC3b(LinkedList<Truple> dataSet) {
		double result = 0;
		for (int i = 0; i < dataSet.size(); i++) {
			feedPattern3b(i, dataSet);
			iterateRight3b();
			result += Math.abs(signum(output3b) - dataSet.get(i).z);
		}
		result /= 2*dataSet.size();
		return result;
	}

	/**
	 * Returns the value of the derivation of our activation function tanh(Beta * b).
	 * @param b3a input parameter.
	 * @return the value of the derivation of our activation function tanh(Beta * b).
	 */
	private double activation3aPrime(double b3a) {
		return BETA3/Math.pow(Math.cosh(BETA3*b3a),2);
	}

	/**
	 * EXERCISE 3
	 * Feed the indicated pattern to the percpetron's input nodes.
	 * @param index The index of the pattern to be fed
	 * @param list The list to select a pattern from
	 */
	private void feedPattern3a(int index, LinkedList<Truple> list) {
		if (index >= list.size()) {
			index = list.size() - 1;
		}
		inputs3a[0] = list.get(index).x;
		inputs3a[1] = list.get(index).y;
	}
	
	/**
	 * EXERCISE 3
	 * Reads the files containing our validation and training sets,
	 * and stores the results in a list.
	 * @param path The path of the file to be read.
	 */
	private LinkedList<Truple> fileToList(String path) throws IOException {
		LinkedList<Truple> list = new LinkedList<>();
		List<String> strList = Files.readAllLines(Paths.get(path), Charset.forName("UTF-8"));
		
		for (String str : strList) {
			//The values are tab separated
			String[] values = str.split("\t");
			list.add(new Truple(
					Double.parseDouble(values[0]),
					Double.parseDouble(values[1]),
					Double.parseDouble(values[2])));
		}
		
		return list;
	}
	
	/**
	 * EXERCISE 3
	 * Performs normalization such that we have zero mean and unit variance
	 * in our training and validation sets
	 * @param list The list to normalize
	 */
	private void normalizeInputs(LinkedList<Truple> list) {
		double meanX = 0;
		double meanY = 0;
		double meanOfSquareX = 0;
		double meanOfSquareY = 0;
		
		//Calculate mean of our two columns
		for (Truple t : list) {
			meanX += t.x;
			meanY += t.y;
		}
		meanX /= list.size();
		meanY /= list.size();
		
		//Subtract the mean from all values in both columns, setting mean to 0
		for (Truple t : list) {
			t.x -= meanX;
			t.y -= meanY;
		}
		
		//Sum up the squares of all our values
		for (Truple t : list) {
			meanOfSquareX += t.x * t.x;
			meanOfSquareY += t.y * t.y;
		}
		
		double varOfX = meanOfSquareX / list.size();
		double varOfY = meanOfSquareY / list.size();
		
		//Divide by the standard deviation to achieve unit variance
		for (Truple t : list) {
			t.x /= Math.sqrt(varOfX);
			t.y /= Math.sqrt(varOfY);
		}
		
	}
	
	/**
	 * EXERCISE 3
	 * Calculates the local field.
	 */
	private double b3a() {
		double result = 0;
		for (int i = 0; i < 2; i++) {
			result += weights3a[i] * inputs3a[i];
		}
		result -= threshold3a;
		return result;
	}
	
	/**
	 * EXERCISE 3
	 * Activation function.
	 */
	private double activation3a(double b) {
		return Math.tanh(BETA3 * b);
	}
	
	public MainNetwork() {
		//Exercise 1
		//task1b();
		
		//Exercise 2
		/*try {
			task2();
		} catch (IOException e) {
			System.out.println("TASK 2A ERROR");
			e.printStackTrace();
		}*/
		
		//Exercise 3
		try {
			//task3a();
			task3b();
		} catch (IOException e) {
			System.out.println("TASK 3 ERROR");
			e.printStackTrace();
		}
	}

	double[] inputs3b = new double[2];
	double[][] hiddenWeights = new double[4][2];
	double[] hiddenlayer = new double[4];
	double[] hiddenThresholds = new double[4];
	double[] outputWeights = new double[4];
	double outputThreshold = 0;
	double output3b = 0;

	private void task3b() throws IOException {
		//Prepare our data
		trainingSet = fileToList("training.txt");
		validationSet = fileToList("validation.txt");
		normalizeInputs(trainingSet);
		normalizeInputs(validationSet);
		double[] classificationsTraining = new double[10];
		double[] classificationsValidation = new double[10];

		//Reset exported data.
		Files.write(Paths.get("3bOutput.txt"), "".getBytes());

		for (int run = 1; run <= 10; run++) {
			trainingEnergyList.clear();
			validationEnergyList.clear();
			System.out.println(".");
			//Set up weights and threshold randomly
			for (int i = 0; i < outputWeights.length; i++) {
				outputWeights[i] = rand.nextDouble() * 0.4 - 0.2;
			}
			for (int i = 0; i < hiddenWeights.length; i++) {
				for (int j = 0; j < hiddenWeights[0].length; j++) {
					hiddenWeights[i][j] = rand.nextDouble() * 0.4 - 0.2;
				}
			}
			outputThreshold = rand.nextDouble() * 2 - 1;
			for (int i = 0; i < hiddenThresholds.length; i++) {
				hiddenThresholds[i] = rand.nextDouble() * 2 - 1;
			}

			//Run the network for a 10^6 iterations
			for (int iteration = 0; iteration < 1000 * 1000; iteration++) {
				int randNum = rand.nextInt(trainingSet.size());
				feedPattern3b(randNum, trainingSet);
				//Iterate Right
				iterateRight3b();

				//Calc errors
				double iError = (trainingSet.get(randNum).z-output3b)*activation3aPrime(bi3b());
				double[] jErrors = new double[4];
				for (int j = 0; j < hiddenlayer.length; j++) {
					jErrors[j] = iError * outputWeights[j] * activation3aPrime(bj3b(j));
				}

				//Update hidden layer weights
				for (int j = 0; j < hiddenWeights.length; j++) {
					for (int k = 0; k < hiddenWeights[0].length; k++) {
						hiddenWeights[j][k] += ETA * jErrors[j] * inputs3b[k];
					}
				}
				//Update hidden thresholds
				for (int j = 0; j < hiddenThresholds.length; j++) {
					hiddenThresholds[j] += -ETA * jErrors[j];
				}
				//Update output weights
				for (int j = 0; j < outputWeights.length; j++) {
					outputWeights[j] += ETA * iError * hiddenlayer[j];
				}
				//Update output threshold
				outputThreshold += -ETA*iError;

				if (iteration % 5000 == 0){
					double trainingEnergy = getEnergy3b(trainingSet);
					trainingEnergyList.add(new Tuple(iteration, trainingEnergy));
					double validationEnergy = getEnergy3b(validationSet);
					validationEnergyList.add(new Tuple(iteration, validationEnergy));
				}
			}
			exportEnergyData(trainingEnergyList, rand.nextDouble() + "," + rand.nextDouble() + ", 0.5", (2*run-1), "3bOutput.txt");
			exportEnergyData(validationEnergyList, "0.5," + rand.nextDouble() + "," + rand.nextDouble(), (2*run),"3bOutput.txt");
			classificationsTraining[run-1] = getC3b(trainingSet);
			classificationsValidation[run-1] = getC3b(validationSet);
		}
		//More output for Mathematica
		StringBuilder sb = new StringBuilder();
		sb.append("Show[");
		for (int i = 1; i <= 10; i++) {
			sb.append("g" + (2*i-1) + ",");
		}
		sb.append(" PlotRange->{{0, 1000000},{0,200}}]\n");
		Files.write(Paths.get("3bOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);
		sb = new StringBuilder();
		sb.append("Show[");
		for (int i = 1; i <= 10; i++) {
			sb.append("g" + (2*i) + ",");
		}
		sb.append(" PlotRange->{{0, 1000000},{0,200}}]\n");
		Files.write(Paths.get("3bOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);

		double averageClass = 0;
		double minimumClass = 400000;
		double varianceClass = 0;
		for (int i = 0; i < 10; i++) {
			averageClass += classificationsTraining[i];
			if (classificationsTraining[i] < minimumClass){
				minimumClass = classificationsTraining[i];
			}
			varianceClass += Math.pow(classificationsTraining[i],2);
		}
		averageClass = (averageClass / 10);
		varianceClass = (varianceClass/10) - Math.pow(averageClass,2);
		System.out.println("Trainingset values:\n Avr:" + averageClass + "\nmin: " + minimumClass + "\nvar:" + varianceClass);

		averageClass = 0;
		minimumClass = 400000;
		varianceClass = 0;
		for (int i = 0; i < 10; i++) {
			averageClass += classificationsValidation[i];
			if (classificationsValidation[i] < minimumClass){
				minimumClass = classificationsValidation[i];
			}
			varianceClass += Math.pow(classificationsValidation[i],2);
		}
		averageClass = (averageClass / 10);
		varianceClass = (varianceClass/10) - Math.pow(averageClass,2);
		System.out.println("Validationset values:\n Avr:" + averageClass + "\nmin: " + minimumClass + "\nvar:" + varianceClass);

		printScatterText(trainingSet, validationSet);
	}

	private void printScatterText(LinkedList<Truple> trainingSet, LinkedList<Truple> validationSet) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("f1 = ListPlot[{");
		for (int i = 0; i < trainingSet.size(); i++) {
			if (trainingSet.get(i).z == 1){
				sb.append("{" + trainingSet.get(i).x + "," + trainingSet.get(i).y + "},");
			}
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("}, PlotStyle->RGBColor[0,0,0]];\n\n");

		sb.append("f2 = ListPlot[{");
		for (int i = 0; i < validationSet.size(); i++) {
			if (validationSet.get(i).z == 1) {
				sb.append("{" + validationSet.get(i).x + "," + validationSet.get(i).y + "},");
			}
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("}, PlotStyle->RGBColor[0.3,0.3,0.3]];\n\n");

		sb.append("f3 = ListPlot[{");
		for (int i = 0; i < trainingSet.size(); i++) {
			if (trainingSet.get(i).z == -1){
				sb.append("{" + trainingSet.get(i).x + "," + trainingSet.get(i).y + "},");
			}
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("}, PlotStyle->RGBColor[0,0,1]];\n\n");

		sb.append("f4 = ListPlot[{");
		for (int i = 0; i < validationSet.size(); i++) {
			if (validationSet.get(i).z == -1) {
				sb.append("{" + validationSet.get(i).x + "," + validationSet.get(i).y + "},");
			}
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("}, PlotStyle->RGBColor[0.5,0.5,1]];\n\n");

		for (int j = 0; j < hiddenlayer.length; j++) {
			sb.append("l"+ j + " = Plot[("+ hiddenWeights[j][0] + "x - " + hiddenThresholds[j] + ")/" + hiddenWeights[j][1] + ", {x,-2,2}];\n");
		}

		sb.append("Show[f2,f1,f3,f4 l1,l2,l3,l0]");
		Files.write(Paths.get("3bOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);

	}

	private double getEnergy3b(LinkedList<Truple> list) {
		double result = 0;
		for (int i = 0; i < list.size(); i++) {
			feedPattern3b(i, list);
			iterateRight3b();
			result += Math.pow(list.get(i).z - output3b, 2);
		}
		result /= 2;
		return result;
	}

	private void iterateRight3b() {
		for (int j = 0; j < hiddenlayer.length; j++) {
			hiddenlayer[j] = activation3a(bj3b(j));
		}
		output3b = activation3a(bi3b());
	}

	private double bi3b() {
		double result = 0;
		for (int j = 0; j < outputWeights.length; j++) {
			result += outputWeights[j]*hiddenlayer[j];
		}
		return result - outputThreshold;
	}

	private double bj3b(int j) {
		double result = 0;
		for (int k = 0; k < inputs3b.length; k++) {
			result += inputs3b[k]*hiddenWeights[j][k];
		}
		return result - hiddenThresholds[j];
	}


	public void feedPattern3b(int index, LinkedList<Truple> list) {
		inputs3b[0] = list.get(index).x;
		inputs3b[1] = list.get(index).y;
	}

	public static void main(String[] args) {
		MainNetwork mn = new MainNetwork();
	}
}
