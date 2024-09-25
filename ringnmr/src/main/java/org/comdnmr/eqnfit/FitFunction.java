/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.comdnmr.eqnfit;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NotPositiveException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.comdnmr.fit.FitQuality;
import org.comdnmr.util.CoMDOptions;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 *
 * @author Bruce Johnson
 */
public abstract class FitFunction implements MultivariateFunction {

    // fixme is there a thread safe RandomGenerator
    public final RandomGenerator DEFAULT_RANDOMGENERATOR = new MersenneTwister(1);
    CoMDOptions options;
    int reportAt = 10;

    EquationType equation;
    long startTime = 0;
    long endTime = 0;
    long fitTime = 0;
    double[][] xValues;
    double[] fieldValues;
    double[] yValues;
    double[] errValues;
    int[] idNums;
    int[][] map;
    int nID = 1;
    boolean reportFitness = false;
    final boolean absMode;
    private static boolean calcError = true;
    double[][] parValues;
    double[] lowerBounds;
    double[] upperBounds;
    final boolean weightFit;

    public class Checker extends SimpleValueChecker {

        public Checker(double relativeThreshold, double absoluteThreshold, int maxIter) {
            super(relativeThreshold, absoluteThreshold, maxIter);
        }

        @Override
        public boolean converged(final int iteration, final PointValuePair previous, final PointValuePair current) {
            boolean converged = super.converged(iteration, previous, current);
            if (reportFitness) {
                if (converged) {
                    System.out.println(previous.getValue() + " " + current.getValue());
                }
                if (converged || (iteration == 1) || ((iteration % reportAt) == 0)) {
                    long time = System.currentTimeMillis();
                    long deltaTime = time - startTime;
                    System.out.println(deltaTime + " " + iteration + " " + current.getValue());
                }
            }
            return converged;
        }
    }

    public FitFunction(CoMDOptions options) {
        this.options = options;
        absMode = options.getAbsValueFit();
        weightFit = options.getWeightFit();
    }

    public abstract void setEquation(String eqName);

    public PointValuePair refine(double[] guess, double[] lowerBounds, double[] upperBounds, double inputSigma, String type) {
        if (type.equals("BOBYQA")) {
            return refineBOBYQA(guess, lowerBounds, upperBounds, inputSigma);
        } else {
            return refineCMAES(guess, lowerBounds, upperBounds, inputSigma);

        }
    }

    public PointValuePair refineCMAES(double[] guess, double[] lowerBounds, double[] upperBounds, double inputSigma) {
        this.lowerBounds = lowerBounds.clone();
        this.upperBounds = upperBounds.clone();
        startTime = System.currentTimeMillis();
        DEFAULT_RANDOMGENERATOR.setSeed(1);
        double lambdaMul = 3.0;
        int lambda = (int) (lambdaMul * FastMath.round(4 + 3 * FastMath.log(guess.length)));
        //int nSteps = guess.length*1000;
        int nSteps = 2000;
        double stopFitness = 0.0;
        int diagOnly = 0;
        double tol = options.getTolerance();
        tol = Math.pow(10.0, tol);
        double[] normLower = new double[guess.length];
        double[] normUpper = new double[guess.length];
        double[] sigma = new double[guess.length];
        Arrays.fill(normLower, 0.0);
        Arrays.fill(normUpper, 100.0);
        Arrays.fill(sigma, inputSigma);
        double[] normGuess = normalize(guess);
        fixGuesses(normGuess);

        //new Checker(100 * Precision.EPSILON, 100 * Precision.SAFE_MIN, nSteps));
        CMAESOptimizer optimizer = new CMAESOptimizer(nSteps, stopFitness, true, diagOnly, 0,
                DEFAULT_RANDOMGENERATOR, true,
                new CPMGFitFunction.Checker(tol, tol, nSteps));
        PointValuePair result = null;

        try {
            result = optimizer.optimize(
                    new CMAESOptimizer.PopulationSize(lambda),
                    new CMAESOptimizer.Sigma(sigma),
                    new MaxEval(2000000),
                    new ObjectiveFunction(this), GoalType.MINIMIZE,
                    new SimpleBounds(normLower, normUpper),
                    new InitialGuess(normGuess));
        } catch (DimensionMismatchException | NotPositiveException | NotStrictlyPositiveException |
                 TooManyEvaluationsException e) {
            e.printStackTrace();
        }
        endTime = System.currentTimeMillis();
        fitTime = endTime - startTime;
        PointValuePair deNormResult = new PointValuePair(deNormalize(result.getPoint()), result.getValue());

        return deNormResult;
    }

    public PointValuePair refineBOBYQA(double[] guess, double[] lowerBounds, double[] upperBounds, double inputSigma) {
        this.lowerBounds = lowerBounds.clone();
        this.upperBounds = upperBounds.clone();
        startTime = System.currentTimeMillis();
        DEFAULT_RANDOMGENERATOR.setSeed(1);
        double lambdaMul = 3.0;
        int lambda = (int) (lambdaMul * FastMath.round(4 + 3 * FastMath.log(guess.length)));
        //int nSteps = guess.length*1000;
        int nSteps = 2000;
        double stopFitness = 0.0;
        int diagOnly = 0;
        double[] normLower = new double[guess.length];
        double[] normUpper = new double[guess.length];
        Arrays.fill(normLower, 0.0);
        Arrays.fill(normUpper, 100.0);
        double[] normGuess = normalize(guess);
        fixGuesses(normGuess);
        //new Checker(100 * Precision.EPSILON, 100 * Precision.SAFE_MIN, nSteps));

        int n = guess.length;
        int nInterp = 2 * n + 1;
        double initialRadius = inputSigma;
        double stopRadius = options.getFinalRadius();
        stopRadius = Math.pow(10.0, stopRadius);

        BOBYQAOptimizer optimizer = new BOBYQAOptimizer(nInterp, initialRadius, stopRadius);
        PointValuePair result = null;

        try {
            result = optimizer.optimize(
                    new MaxEval(2000000),
                    new ObjectiveFunction(this), GoalType.MINIMIZE,
                    new SimpleBounds(normLower, normUpper),
                    new InitialGuess(normGuess));
        } catch (Exception e) {
            e.printStackTrace();
        }
        endTime = System.currentTimeMillis();
        fitTime = endTime - startTime;
        PointValuePair deNormResult = new PointValuePair(deNormalize(result.getPoint()), result.getValue());

        return deNormResult;
    }

    void fixGuesses(double[] guesses) {
        for (int i = 0; i < guesses.length; i++) {
            if (guesses[i] > 98.0) {
                guesses[i] = 98.0;
            } else if (guesses[i] < 2) {
                guesses[i] = 2.0;
            }
        }
    }

    double[] normalize(double[] pars) {
        double[] normPars = new double[pars.length];
        for (int i = 0; i < pars.length; i++) {
            normPars[i] = 100.0 * (pars[i] - lowerBounds[i]) / (upperBounds[i] - lowerBounds[i]);
        }
        return normPars;
    }

    double[] deNormalize(double[] normPars) {
        double[] pars = new double[normPars.length];
        for (int i = 0; i < pars.length; i++) {
            pars[i] = normPars[i] / 100.0 * (upperBounds[i] - lowerBounds[i]) + lowerBounds[i];
        }
        return pars;
    }

    public double[] simY(double[] par) {
        double[] yCalc = equation.calculate(par, map[0], xValues, idNums[0]);
        return yCalc;
    }

    public static void setCalcError(boolean state) {
        calcError = state;
    }

    public static boolean getCalcError() {
        return calcError;
    }

    public abstract int[] getMask();

    public abstract double[] simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options);

    public abstract double[] simBoundsStream(double[] start,
                                             double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options);

    public abstract double[][] getSimPars();

    public void setXY(double[][] x, double[] y) {
        this.xValues = x;
        this.yValues = y;
        this.idNums = new int[x[0].length];
    }

    public void setErr(double[] err) {
        this.errValues = err;
    }

    public void setFieldValues(double[] fieldValues) {
        this.fieldValues = fieldValues;
    }

    public double[] getFieldValues() {
        return fieldValues;
    }

    final boolean setNID() {
        nID = Arrays.stream(idNums).max().getAsInt() + 1;
        boolean[] checkIDs = new boolean[nID];
        Arrays.stream(idNums).forEach(id -> checkIDs[id] = true);
        return IntStream.range(0, checkIDs.length).anyMatch(id -> checkIDs[id] == false);
    }

    public void setIds(int[] idNums) throws IllegalArgumentException {
        this.idNums = idNums;
        if (setNID()) {
            //for (int id : idNums) {
            //System.out.print(id + " ");
            //}
            //System.out.println("");
            throw new IllegalArgumentException("Invalid idNums, some values not used in setIds");
        }
    }

    public double[] guess() {
        double[] guess = equation.guess(xValues, yValues, map, idNums, nID);
        return guess;
    }

    public double[][] boundaries(double[] guesses) {
        //return equation.boundaries(xValues, yValues, fieldValues[0]);
        double[][] boundaries = equation.boundaries(guesses, xValues, yValues, map, idNums, nID);
        return boundaries;
    }

    public double calculate(double[] par, int[] map, double[] x, int idNum) {
        return equation.calculate(par, map, x, idNum);
    }

    public double getRSS(double[] par) {
        double rss = 0.0;
        double[] x = new double[xValues.length];
        for (int i = 0; i < yValues.length; i++) {
            for (int j = 0; j < x.length; j++) {
                x[j] = xValues[j][i];
            }
            final double value;
            value = calculate(par, map[idNums[i]], x, idNums[i]);
            double delta = value - yValues[i];
            rss += delta * delta;
        }
        return rss;
    }

    public double getRMS(double[] par) {
        double rss = 0.0;
        double[] x = new double[xValues.length];
        for (int i = 0; i < yValues.length; i++) {
            for (int j = 0; j < x.length; j++) {
                x[j] = xValues[j][i];
            }
            final double value;
            value = calculate(par, map[idNums[i]], x, idNums[i]);

            double delta = value - yValues[i];
            rss += delta * delta;
        }
        return Math.sqrt(rss / yValues.length);
    }

    public double getAIC(double[] par) {
        double rss = getRSS(par);
        int k = par.length;
        int n = yValues.length;
        return 2 * k + n * Math.log(rss);
    }

    public double getAICc(double[] par) {
        double rss = getRSS(par);
        int k = par.length;
        int n = yValues.length;
        double aic = getAIC(par);
        double aicc = aic + 2 * k * (k + 1) / (n - k - 1);
        return aicc;
    }

    public FitQuality getFitQuality(double[] pars) {
        return new FitQuality(getRMS(pars), getAIC(pars),
                getAICc(pars), getReducedChiSq(pars), pars.length);
    }

    public double getReducedChiSq(double[] par) {
        double rss = 0.0;
        double[] x = new double[xValues.length];
        for (int i = 0; i < yValues.length; i++) {
            for (int j = 0; j < x.length; j++) {
                x[j] = xValues[j][i];
            }
            final double value;
            value = calculate(par, map[idNums[i]], x, idNums[i]);

            double delta = (value - yValues[i]) / errValues[i];
            rss += delta * delta;
        }
        return rss / (yValues.length - par.length);
    }

    public double[] getPredicted(double[] par) {
        double[] yPred = simY(par);
        return yPred;
    }

    public void setMap(int[][] map) {
        this.map = map;
    }

    public void setMap(int[] stateCount, int[][] states) {
//        System.out.println("states ");
        //for (int i = 0; i < states.length; i++) {
            //for (int j = 0; j < states[i].length; j++) {
                //System.out.print(states[i][j] + " ");
            //}
            //System.out.println("");
        //}
//        System.out.println(" stateCount");
//        for (int j = 0; j < stateCount.length; j++) {
//            System.out.print(stateCount[j] + " ");
//        }
//        System.out.println(" ");
        this.map = new int[states.length][stateCount.length];
        this.map = equation.makeMap(stateCount, states, getMask());
//        System.out.println("map ");
//        for (int i = 0; i < map.length; i++) {
//            for (int j = 0; j < map[i].length; j++) {
//                System.out.print(map[i][j] + " ");
//            }
//            System.out.println("");
//        }
//        System.out.println("done ");
    }

    public int[][] getMap() {
        return map;
    }

    public String[] getParNames() {
        return equation.getParNames();
    }

    public int getNGroupPars() {
        return equation.getNGroupPars();
    }

    public boolean checkID(int[] idValues, int[] testValues, int checkNum) {
        int n = Arrays.stream(idValues).max().getAsInt() + 1;
        int[] count = new int[n];
        for (int i = 0; i < testValues.length; i++) {
            count[testValues[i]]++;
        }
        boolean ok = true;
        for (int i = 0; i < count.length; i++) {
            if (count[i] < checkNum) {
                ok = false;
                break;
            }
        }
        return ok;
    }

/*    public void dump(double[] par, double[][] xValues, double field) {
        double[] x = new double[xValues.length];
        for (int i = 0; i < yValues.length; i++) {
            for (int j = 0; j < x.length; j++) {
                x[j] = xValues[j][i];
            }
            double yCalc = equation.calculate(par, map[idNums[i]], x, idNums[i]);
            System.out.printf("%8.5f %8.5f %8.5f\n", x[0], yCalc, field);
        }
    }

    public void dump(double[] par) {
        for (int i = 0; i < xValues.length; i++) {
            double yCalc = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], fieldValues[i]);
            System.out.printf("%8.5f %8.5f %8.5f\n", xValues[0][i], yValues[i], yCalc);

        }

        DescriptiveStatistics dstat = new DescriptiveStatistics(xValues[0]);
        double min = dstat.getMin();
        double max = dstat.getMax();
        int n = 20;
        double[] ax = new double[1];
        double delta = (max - min) / (n - 1);
        for (int i = 0; i < n; i++) {
            double x = min + i * delta;
            ax[0] = x;
            double y = equation.calculate(par, map[idNums[i]], ax, idNums[i], fieldValues[i]);
            System.out.printf("%8.5f %8.5f\n", x, y);
        }
    }

    void dumpValues() {
        for (int i = 0; i < yValues.length; i++) {
            System.out.printf("%3d %8.3f %8.3f %8.3f %8.3f %3d\n", i, xValues[0][i], yValues[i], errValues[i], fieldValues[i], idNums[i]);
        }
    }*/
}
