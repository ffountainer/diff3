// packages for the graph
package org.example;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import javax.swing.*;


import java.sql.SQLOutput;
import java.util.ArrayList;


public class Main {
    public static void main(String[] args) {
        solve("Task", 1, 1, 0);
    }

    public static void solve(String task, double x_0, double y_0, double t_0_) {
        ArrayList<ArrayList<Pair>> secondOrderAnswers = new ArrayList<>();
        ArrayList<ArrayList<Double>> time = new ArrayList<>();
        double[] deltaT = {0.00001, 0.000005};
        for (int delta_index = 0; delta_index < deltaT.length; delta_index++) {
            secondOrderAnswers.add(new ArrayList<>());

            double x = x_0;
            double y = y_0;
            double t = t_0_;

            ArrayList<Double> timeCur = new ArrayList<>();

            while (t <= 10.01) {
                Triple triple = secondOrderRungeKuttaIteration(y, x, t, delta_index, 1, deltaT);
                t = triple.firstElement;
                x = triple.secondElement;
                y = triple.thirdElement;
                Pair pair = new Pair(x, y);
                timeCur.add(t);
                secondOrderAnswers.get(delta_index).add(pair);

            }

            time.add(timeCur);

        }


        ArrayList<ArrayList<Double>> secondOrderErrors = new ArrayList<>();


        for (int i = 1; i < secondOrderAnswers.size(); i++) {
            secondOrderErrors.add(new ArrayList<>());

            int step2 = (int) ((int)Math.round(10/deltaT[i]));
            int step1 = (int) ((int)Math.round(10/deltaT[i - 1]));

            int ratio2 = step2 / 10;
            int ratio1 = step1 / 10;

            int k = ratio2;
            int j = ratio1;


            for (int counter = 0; counter < 10; counter++) {

                double elementx1 = secondOrderAnswers.get(i - 1).get(j).getFirstElement();
                double elementx2 = secondOrderAnswers.get(i).get(k).getFirstElement();


                double elementy1 = secondOrderAnswers.get(i - 1).get(j).getSecondElement();
                double elementy2 = secondOrderAnswers.get(i).get(k).getSecondElement();


                secondOrderErrors.get(i - 1).add(Math.abs(elementx2 - elementx1) + Math.abs(elementy2 - elementy1));

                j += ratio1;
                k += ratio2;
            }
        }

        int ind = 0;
        int num = -1;
        for (ArrayList<Double> errors : secondOrderErrors) {
            double absError = errors.get(0);
            if (absError < 1e-8) {
                num = ind;
            }
            for (double error : errors) {
                if (error > absError && error <  1e-8) {
                    num = ind;
                }
                absError = Math.max(absError, error);
            }

            ind += 1;
        }

        int[] N = new int[deltaT.length];
        for (int i = 0; i< deltaT.length; i++) {
            N[i] = (int) ((int)Math.round(10/deltaT[i]));
        }

        System.out.println();
        System.out.println("Second order");

        if (num == -1) {
            num = 0;
            System.out.println("The max error is more than 10^-8");
        }

        System.out.println("Local errors for 2nd order Runge-Kutta:");
        System.out.println("N = " + N[num]);
        for (double error : secondOrderErrors.get(num)) {
            System.out.print(error + " ");
        }

        System.out.println();

        System.out.println();


        ArrayList<ArrayList<Pair>> fourthOrderAnswers = new ArrayList<>();

        for (int delta_index = 0; delta_index < deltaT.length; delta_index++) {
            fourthOrderAnswers.add(new ArrayList<>());

            double x = x_0;
            double y = y_0;
            double t = t_0_;

            ArrayList<Double> timeCur = new ArrayList<>();

            while (t <= 10.01) {
                Triple triple = fourthOrderRungeKuttaIteration(y, x, t, delta_index, 1, deltaT);
                t = triple.firstElement;
                x = triple.secondElement;
                y = triple.thirdElement;
                Pair pair = new Pair(x, y);
                timeCur.add(t);
                fourthOrderAnswers.get(delta_index).add(pair);

            }

            time.add(timeCur);

        }



        ArrayList<ArrayList<Double>> fourthOrderErrors = new ArrayList<>();

        for (int i = 1; i < fourthOrderAnswers.size(); i++) {
            fourthOrderErrors.add(new ArrayList<>());

            int step2 = (int) ((int)Math.round(10/deltaT[i]));
            int step1 = (int) ((int)Math.round(10/deltaT[i - 1]));

            int ratio2 = step2 / 10;
            int ratio1 = step1 / 10;

            int k = ratio2;
            int j = ratio1;


            for (int counter = 0; counter < 10; counter++) {

                double elementx1 = fourthOrderAnswers.get(i - 1).get(j).getFirstElement();
                double elementx2 = fourthOrderAnswers.get(i).get(k).getFirstElement();


                double elementy1 = fourthOrderAnswers.get(i - 1).get(j).getSecondElement();
                double elementy2 = fourthOrderAnswers.get(i).get(k).getSecondElement();


                fourthOrderErrors.get(i - 1).add(Math.abs(elementx2 - elementx1) + Math.abs(elementy2 - elementy1));

                j += ratio1;
                k += ratio2;
            }
        }

        int ind4 = 0;
        int num4 = -1;
        for (ArrayList<Double> errors : fourthOrderErrors) {
            double absError = errors.get(0);
            if (absError < 1e-8) {
                num4 = ind4;
            }
            for (double error : errors) {
                if (error > absError && error <  1e-8) {
                    num4 = ind4;
                }
                absError = Math.max(absError, error);
            }

            ind += 1;
        }

        int[] N4 = new int[deltaT.length];
        for (int i = 0; i< deltaT.length; i++) {
            N4[i] = (int) ((int)Math.round(10/deltaT[i]));
        }


        System.out.println("Fourth order");

        if (num4 == -1) {
            num4 = 0;
            System.out.println("The max error is more than 10^-8");
        }

        System.out.println("Local errors for 4nd order Runge-Kutta:");
        System.out.println("N = " + N[num4]);
        for (double error : fourthOrderErrors.get(num4)) {
            System.out.print(error + " ");
        }


        System.out.println();



        double[] secondOrderDataY = new double[secondOrderAnswers.get(num).size()];
        for (int i = 0; i < secondOrderAnswers.get(num).size(); i++) {
            secondOrderDataY[i] = secondOrderAnswers.get(num).get(i).getSecondElement();
        }


        double[] secondOrderDataX = new double[secondOrderAnswers.get(num).size()];
        for (int i = 0; i < secondOrderAnswers.get(num).size(); i++) {
            secondOrderDataX[i] = secondOrderAnswers.get(num).get(i).getFirstElement();
        }


        XYChart chart = new XYChartBuilder()
                .title("VELE67I03")
                .xAxisTitle("x(t)")
                .yAxisTitle("y(t)")
                .build();


        chart.addSeries("2nd order Runge-Kutta", secondOrderDataX, secondOrderDataY);

        JPanel chartPanel = new XChartPanel<>(chart);
        JFrame frame = new JFrame(task + " - Veronika Levasheva");
        frame.add(chartPanel);
        frame.setSize(800, 600);
        frame.show();

        double[] fourthOrderDataY = new double[fourthOrderAnswers.get(num4).size()];
        for (int i = 0; i < fourthOrderAnswers.get(num4).size(); i++) {
            fourthOrderDataY[i] = fourthOrderAnswers.get(num4).get(i).getSecondElement();
        }

        double[] fourthOrderDataX = new double[fourthOrderAnswers.get(num4).size()];
        for (int i = 0; i < fourthOrderAnswers.get(num4).size(); i++) {
            fourthOrderDataX[i] = fourthOrderAnswers.get(num4).get(i).getFirstElement();
        }

        XYChart chart2 = new XYChartBuilder()
                .title("VELE67I03")
                .xAxisTitle("x(t)")
                .yAxisTitle("y(t)")
                .build();


        chart2.addSeries("4nd order Runge-Kutta", fourthOrderDataX, fourthOrderDataY);

        JPanel chartPanel2 = new XChartPanel<>(chart2);
        JFrame frame2 = new JFrame(task + " - Veronika Levasheva");
        frame2.add(chartPanel2);
        frame2.setSize(800, 600);
        frame2.show();
    }

    public static Triple secondOrderRungeKuttaIteration(double y_t, double x_t, double t, int deltaTInd, int direction, double[] deltaT) {

        double h = deltaT[deltaTInd];
//        System.out.println("h " + h + " t " + t + " x " + x_t + " y "  + y_t);

        double dy_1 = 1.0;
        double dy_2 = 1.0;

        double dx_1 = 1.0;
        double dx_2 = 1.0;


        dx_1 = h * (-(x_t*x_t) - y_t);
        dy_1 = h * (2*x_t - y_t);

        dx_2 = h * (- (x_t + dx_1)*(x_t + dx_1) - (y_t + dy_1));
        dy_2 = h * ( 2 * (x_t + dx_1) - (y_t + dy_1) );

        x_t = x_t + (dx_1 + dx_2) / 2;
        y_t = y_t + (dy_1 + dy_2) / 2;

        t = t + h;

        return new Triple(t, x_t, y_t);

    }

    public static Triple fourthOrderRungeKuttaIteration(double y_t, double x_t, double t, int deltaTInd, int direction, double[] deltaT) {
        deltaT[deltaTInd] = deltaT[deltaTInd] * direction;

        double dy_1 = 1.0;
        double dy_2 = 1.0;
        double dy_3 = 1.0;
        double dy_4 = 1.0;

        double dx_1 = 1.0;
        double dx_2 = 1.0;
        double dx_3 = 1.0;
        double dx_4 = 1.0;

        double h = deltaT[deltaTInd];

        dx_1 = -(x_t*x_t) - y_t;
        dy_1 = 2*x_t - y_t;

        dx_2 = - (x_t + (h/2) * dx_1)*(x_t + (h/2) * dx_1) - (y_t + (h/2) * dy_1);
        dy_2 = 2 * (x_t + (h/2) * dx_1) - (y_t + (h/2) * dy_1);

        dx_3 = - (x_t + (h/2) * dx_2)*(x_t + (h/2) * dx_2) - (y_t + (h/2) * dy_2);
        dy_3 = 2 * (x_t + (h/2) * dx_2) - (y_t + (h/2) * dy_2);

        dx_4 = - (x_t + h * dx_3)*(x_t + h * dx_3) - (y_t + h * dy_3);
        dy_4 = 2 * (x_t + h * dx_3) - (y_t + h * dy_3);

        x_t = x_t + (h/6) * (dx_1 + 2 * dx_2 + 2 * dx_3 + dx_4);
        y_t = y_t + (h/6) * (dy_1 + 2 * dy_2 + 2 * dy_3 + dy_4);

        t = t + h;


        return new Triple(t, x_t, y_t);

    }
}

class Triple{
    double firstElement;
    double secondElement;
    double thirdElement;

    public Triple(double first, double second, double third) {
        this.firstElement = first;
        this.secondElement = second;
        this.thirdElement = third;

    }

    public double getFirstElement() {
        return firstElement;
    }

    public double getSecondElement() {
        return secondElement;
    }

    public double getThirdElement() {
        return thirdElement;
    }
}

class Pair{
    double firstElement;
    double secondElement;

    public Pair(double first, double second) {
        this.firstElement = first;
        this.secondElement = second;

    }

    public double getFirstElement() {
        return firstElement;
    }

    public double getSecondElement() {
        return secondElement;
    }

}
