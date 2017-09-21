import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

public class MainNetwork {
    private final int N = 200;
    private int P = 40;
    private double BETA2 = 2;
    private double ETA = 0.02;
    private double BETA3 = 0.5;
    private double[][] w = new double[N][N];
    private int[][] patterns = new int[P][N];
    private int[] network = new int[N];
    private Random rand = new Random();

    private class Tuple{
        public double x;
        public double y;
        public Tuple(double x, double y){
            this.x = x;
            this.y = y;
        }
    }

    public int signum(double x){
        if (x<0){
            return -1;
        }else{
            return +1;
        }
    }

    public void generatePatterns(){
        for (int i = 0; i < patterns.length; i++) {
            for (int j = 0; j < patterns[0].length; j++) {
                patterns[i][j] = rand.nextBoolean()?1:-1;
            }
        }
    }

    public void calculateWeights(){
        for (int i = 0; i < w.length; i++) {
            for (int j = 0; j < w.length; j++) {
                double sum = 0;
                for (int k = 0; k < P; k++) {
                    sum += patterns[k][i] * patterns[k][j];
                }
                sum = sum/N;
                w[i][j] = sum;
            }
        }
    }

    public void feed(int index){
        //arrays indexed by 0;
        network = new int[N];
        index--;
        for (int i = 0; i < N; i++) {
            network[i] = patterns[index][i];
        }
    }

    private void synchronousDeterministicUpdate(int steps){
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

    private void printPatterns(){
        for (int i = 0; i < P; i++) {
            for (int j = 0; j < N; j++) {
                if (patterns[i][j] == 1){
                    System.out.print("[+]");
                }else{
                    System.out.print("[-]");
                }
            }
            System.out.println("");
        }
    }

    private void printNetwork(){
        System.out.println("Network state:");
        for (int i = 0; i < N; i++) {
            if (network[i] == 1){
                System.out.print("[+]");
            }else{
                System.out.print("[-]");
            }
        }
        System.out.println("");
    }

    public void asynchronousStochasticUpdate(int steps){
        for (int step = 0; step < steps; step++) {
            //Choose a random bit to update
            int i = (int)(rand.nextDouble() * (N));
            if (rand.nextDouble() < g(i)){
                network[i] = 1;
            } else {
                network[i] = -1;
            }
        }
    }

    public double g(int i){
        return 1 / (1 + Math.exp(-2* BETA2 *b(i)));
    }

    public double b(int i){
        double sum = 0;
        for (int j = 0; j < N; j++) {
            sum += w[i][j]*network[j];
        }
        return sum;
    }

    public double m(int pattern){
        double sum = 0;
        for (int i = 0; i < N; i++) {
            sum += patterns[pattern-1][i] * network[i];
        }
        sum = sum / N;
        return sum;
    }

    public MainNetwork() throws IOException {

        //task1b();

        //task2a();

        task3a();
    }

    public void task1b(){
        int faultyUpdates = 0;
        for (int i=0; i<500; i++) {
            //Randomly generated 50/50
            generatePatterns();

            //Hebb's Rule
            calculateWeights();

            //Feed a given pattern.
            feed(1);

            synchronousDeterministicUpdate(1);

            for (int j = 0; j < N; j++) {
                if (patterns[0][j] != network[j]){
                    faultyUpdates++;
                }
            }
            if (i%25 == 0){
                System.out.print(".");
            }
        }
        System.out.println((double)faultyUpdates/100000);
    }

    public void task2a() throws IOException {
        StringBuilder sb;
        Files.write(Paths.get("2aOutput.txt"), "".getBytes());
        int iterations = 20;
        LinkedList<Tuple> steadyList = new LinkedList<>();
        boolean steadyFound = false;
        double maxX = 0;
        for (int i=1;i<=iterations;i++) {
            LinkedList<Tuple> dataList = new LinkedList<>();
            sb = new StringBuilder();
            generatePatterns();
            calculateWeights();
            feed(1);
            int step = -1;
            double m = 0;
            double averageM = 0;
            int BATCH_SIZE = 50;

            steadyFound = false;
            while (true){
                step++;
                asynchronousStochasticUpdate(1);
                m += getM(1);
                if (step%BATCH_SIZE == 0){
                    if (step != 0){
                        m = m / BATCH_SIZE;
                    }
                    averageM = ((averageM * dataList.size())+(m))/(dataList.size()+1);
                    dataList.add(new Tuple(step, averageM));

                    if (dataList.size() > 7 && !steadyFound) {
                        if (Math.abs(dataList.get(dataList.size()-1).y - (dataList.get(dataList.size()-2).y + dataList.get(dataList.size()-3).y + dataList.get(dataList.size()-4).y + dataList.get(dataList.size()-5).y + dataList.get(dataList.size()-6).y + dataList.get(dataList.size()-7).y)/6) < 0.0001){
                            steadyList.add(new Tuple((double)step, 0));
                            steadyFound = true;
                        }
                    }
                    if (step == 100000){
                        maxX = 100000;
                        break;
                    }
                    m=0;
                }
            }
            sb.append("g" + i + " = ListPlot[{");
            for (int j = 0; j < dataList.size(); j++) {
                sb.append("{"+dataList.get(j).x + "," +dataList.get(j).y + "}");
                if (j != dataList.size()-1){
                    sb.append(", ");
                }
            }
            sb.append("}, Joined -> True, PlotRange->{0,1}, PlotStyle->RGBColor["+rand.nextDouble()+","+rand.nextDouble()+","+rand.nextDouble()+"]];\n\n");
            try {
                Files.write(Paths.get("2aOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        double steadySteps = 0;
        for (Tuple t :
                steadyList) {
            System.out.println(t.x);
            steadySteps += t.x;
        }
        steadySteps = steadySteps / steadyList.size();
        System.out.println("Average steady state steps: "+steadySteps);

        sb = new StringBuilder();
        sb.append("Show[");
        for (int i = 1; i <= iterations; i++) {
            sb.append("g" + i + ",");
        }
        sb.append(" PlotRange->{{0, "+maxX+"},{0,1}}]");
        try {
            Files.write(Paths.get("2aOutput.txt"), sb.toString().getBytes(), StandardOpenOption.APPEND);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private double[] weights3a = new double[2];
    private double threshold3a = 0;
    private double[] inputs3a = new double[2];
    private LinkedList<Truple> trainingSet;
    private LinkedList<Truple> validationSet;

    private void task3a() {
        trainingSet = fileToList("training.txt");
        validationSet = fileToList("validation.txt");
        normalizeInputs(trainingSet);
        normalizeInputs(validationSet);

        for (int run = 1; run <= 10; run++) {
            double output = 0;
            for (int i = 0; i < 2; i++) {
                weights3a[i] = rand.nextDouble() * 0.4 - 0.2;
            }
            threshold3a = rand.nextDouble() * 2 - 1;
            for (int i = 0; i < 1000*1000; i++) {
                int randNum = rand.nextInt(trainingSet.size());
                feedPattern3a(randNum, trainingSet);
                output = activation3a(b3a());
                //Update weights
                weights3a[0] += 0;
            }
        }
    }

    private void feedPattern3a(int index, LinkedList<Truple> list) {
        inputs3a[0] = list.get(index).x;
        inputs3a[1] = list.get(index).y;
    }

    private LinkedList<Truple> fileToList(String path) {
        LinkedList<Truple> list = new LinkedList<>();
        try {
            List<String> strList = Files.readAllLines(Paths.get(path), Charset.forName("UTF-8"));
            for (String str :
                    strList) {
                String[] values = str.split("\t");
                list.add(new Truple(Double.parseDouble(values[0]), Double.parseDouble(values[1]), Double.parseDouble(values[2])));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return list;
    }

    private void normalizeInputs(LinkedList<Truple> list){
        double meanX = 0;
        double meanY = 0;
        double meanOfSquareX = 0;
        double meanOfSquareY = 0;

        //Calculate mean of our two columns
        for (Truple trup : list) {
            meanX += trup.x;
            meanY += trup.y;
        }
        meanX /= list.size();
        meanY /= list.size();

        //Subtract the mean from all values in both columns, setting mean to 0
        for (Truple t : list) {
            t.x -= meanX;
            t.y -= meanY;
        }

        //Sum up the squares of all our values.
        for (Truple trup : list) {
            meanOfSquareX += trup.x*trup.x;
            meanOfSquareY += trup.y*trup.y;
        }
        double varOfX = meanOfSquareX / list.size();
        double varOfY = meanOfSquareY / list.size();
        System.out.println("XVAR "+ varOfX);
        //Divide by the standard deviation to achieve unit variance.
        for (Truple t : list) {
            t.x /= Math.sqrt(varOfX);
            t.y /= Math.sqrt(varOfY);
        }

    }

    private double activation3a(double b){
        return Math.tanh(BETA3*b);
    }

    private double b3a(){
        double result = 0;
        for (int i = 0; i < 2; i++) {
            result += weights3a[i]*inputs3a[i];
        }
        result -= threshold3a;
        return result;
    }

    private double getM(int pattern) {
        double result = 0;
        for (int i = 0; i < N; i++) {
            result += network[i]*patterns[pattern-1][i];
        }
        result = result / N;
        return result;
    }

    public static void main(String[] args) throws IOException {
        MainNetwork mn = new MainNetwork();
    }

    private class Truple {
        public double x;
        public double y;
        public double z;
        public Truple(double x, double y, double z){
            this.x = x;
            this.y = y;
            this.z = z;
        }
    }
}
