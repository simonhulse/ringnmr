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

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;
import org.comdnmr.util.CoMDOptions;

public class CPMGFitFunction extends FitFunction {

    static RandomGenerator random = new SynchronizedRandomGenerator(new Well19937c());
    int[] r2Mask = {0, 1, 3};
    double[] rexErrors = new double[nID];

    public CPMGFitFunction(CoMDOptions options) {
        super(options);
        this.equation = CPMGEquation.CPMGFAST;
    }

    @Override
    public void setEquation(String eqName) {
        equation = CPMGEquation.valueOf(eqName.toUpperCase());
    }

    public CPMGFitFunction(CoMDOptions options, double[][] x, double[] y, double[] err) throws IllegalArgumentException {
        this(options, x, y, err, new int[x.length]);
    }

    public CPMGFitFunction(CoMDOptions options, double[][] x, double[] y, double[] err, int[] idNums) throws IllegalArgumentException {
        super(options);
        this.xValues = new double[1][];
        this.xValues[0] = x[0].clone();
        this.yValues = y.clone();
        this.errValues = err.clone();
        this.idNums = idNums.clone();
        this.idNums = new int[yValues.length];
        this.equation = CPMGEquation.CPMGFAST;
        if (setNID()) {
            throw new IllegalArgumentException("Invalid idNums, some values not used");
        }
    }

    public static int getNPars(int[][] map) {
        int maxIndex = 0;
        for (int[] map1 : map) {
            for (int i : map1) {
                maxIndex = Math.max(i, maxIndex);
            }
        }
        return maxIndex + 1;
    }

    @Override
    public int[] getMask() {
        return r2Mask;
    }

    @Override
    public double[][] getSimPars() {
        return parValues;
    }

    @Override
    public double value(double[] normPar) {
        double[] par = deNormalize(normPar);
        double sumAbs = 0.0;
        double sumSq = 0.0;
        double[] ax = new double[4];
        for (int i = 0; i < yValues.length; i++) {
            final double value;
            for (int j = 0; j < xValues.length; j++) {
                ax[j] = xValues[j][i];
            }
            value = equation.calculate(par, map[idNums[i]], ax, idNums[i]);
            double delta = (value - yValues[i]);
            if (weightFit) {
                delta /= errValues[i];
            }
            sumAbs += FastMath.abs(delta);
            sumSq += delta * delta;
        }
        if (reportFitness) {
            double rms = Math.sqrt(sumSq / yValues.length);
            for (double p : par) {
                System.out.printf("%7.3f ", p);
            }
            System.out.printf("%7.3f %7.3f %7.3f\n",
                    sumSq / (yValues.length - par.length), sumAbs / (yValues.length - par.length), rms);
        }

        if (absMode) {
            return sumAbs / (yValues.length - par.length);
        } else {
            return sumSq / (yValues.length - par.length);
        }
    }

    @Override
    public double[] simY(double[] par) {
        double[] ax = new double[4];
        double[] yv = new double[yValues.length];
        for (int i = 0; i < yValues.length; i++) {
            for (int j = 0; j < xValues.length; j++) {
                ax[j] = xValues[j][i];
            }
            yv[i] = equation.calculate(par, map[idNums[i]], ax, idNums[i]);
        }
        return yv;
    }

    public List<Double> simY(double[] par, double[][] xValues) {
        ArrayList<Double> result = new ArrayList<>();
        double[] ax = new double[4];
        for (int i = 0; i < yValues.length; i++) {
            ax[0] = xValues[0][i];
            for (int j = 0; j < xValues.length; j++) {
                ax[j] = xValues[j][i];
            }

            double yCalc = equation.calculate(par, map[idNums[i]], ax, idNums[i]);
            result.add(yCalc);
        }
        return result;
    }

    public static double rDisp(double[] par, double x) {
        double a = par[0];
        double b = par[1];
        double c = par[2];
        // fixme check for x=0;
        return Math.exp(-(a - a * Math.sin(b / x) / (b / x) + c));
    }

    public double[] getRex(double[] pars, double field) {
        double[] result = new double[nID];
        for (int i = 0; i < map.length; i++) {
            result[i] = equation.getRex(pars, map[i],field);
        }
        return result;
    }

    public double getRex(double[] pars, int id, double field) {
        return equation.getRex(pars, map[id], field);
    }

    public double[] getRexError() {
        return rexErrors.clone();
    }

    public double getKex(double[] pars) {
        return equation.getKex(pars);
    }

    @Override
    public Optional<double[]> simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        reportFitness = false;
        int nSim = options.getSampleSize();
        int nPar = start.length;
        parValues = new double[nPar + 1][nSim];
        double[] yPred = simY(start);
        double[] yValuesOrig = yValues.clone();
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        String optimizer = options.getBootStrapOptimizer();

        for (int i = 0; i < nSim; i++) {
            for (int k = 0; k < yValues.length; k++) {
                yValues[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            var resultOpt = refine(start, lowerBounds, upperBounds,
                    inputSigma, optimizer);
            PointValuePair result;
            if (resultOpt.isEmpty()) {
                return Optional.empty();
            } else {
                result = resultOpt.get();
            }

            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
            parValues[nPar][i] = result.getValue();
            if (equation == CPMGEquation.CPMGSLOW) {
                for (int j = 0; j < map.length; j++) {
                    rexValues[j][i] = equation.getRex(result.getPoint(), map[j], fieldValues[0]);
                }
            }
        }
        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        if (equation == CPMGEquation.CPMGSLOW) {
            for (int j = 0; j < nID; j++) {
                DescriptiveStatistics dStat = new DescriptiveStatistics(rexValues[j]);
                rexErrors[j] = dStat.getStandardDeviation();
            }
        }
        yValues = yValuesOrig;
        return Optional.of(parSDev);
    }

    private synchronized CPMGFitFunction setupParametricBootstrap(double[] yPred) {
        double[] newY = new double[yValues.length];
        for (int k = 0; k < yValues.length; k++) {
            newY[k] = yPred[k] + errValues[k] * random.nextGaussian();
        }
        CPMGFitFunction rDisp = new CPMGFitFunction(options, xValues, newY, errValues, idNums);
        rDisp.setEquation(equation.getName());
        rDisp.setXY(xValues, newY);
        rDisp.setIds(idNums);
        rDisp.setMap(map);
        return rDisp;
    }

    private CPMGFitFunction setupNonParametricBootstrap(double[] yPred) {
        CPMGFitFunction rDisp = new CPMGFitFunction(options, xValues, yValues, errValues, idNums);
        rDisp.setEquation(equation.getName());
        double[][] newX = new double[xValues.length][yValues.length];
        double[] newY = new double[yValues.length];
        double[] newErr = new double[yValues.length];
        int[] newID = new int[yValues.length];
        int iTry = 0;
        do {
            for (int k = 0; k < yValues.length; k++) {
                int rI = random.nextInt(yValues.length);
                for (int j=0;j<xValues.length;j++) {
                    newX[j][k] = xValues[j][rI];
                }
                newY[k] = yValues[rI];
                newErr[k] = errValues[rI];
                newID[k] = idNums[rI];
            }
            iTry++;
        } while (!checkID(idNums, newID, 2) && (iTry < 10));
        // fixme  idNum should be set in above loop
        rDisp.setXY(newX, newY);
        rDisp.setErr(newErr);
        rDisp.setIds(newID);
        rDisp.setMap(map);
        return rDisp;
    }

    @Override
    public Optional<double[]> simBoundsStream(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = options.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        double[] yPred = simY(start);
        String optimizer = options.getBootStrapOptimizer();
        AtomicBoolean hadError = new AtomicBoolean(false);
        IntStream.range(0, nSim).parallel().forEach(i -> {
            CPMGFitFunction rDisp;
            if (options.getNonParametricBootstrap()) {
                rDisp = setupNonParametricBootstrap(yPred);
            } else {
                rDisp = setupParametricBootstrap(yPred);
            }

            var resultOpt = rDisp.refine(start, lowerBounds, upperBounds,
                    inputSigma, optimizer);
            if (resultOpt.isEmpty()) {
                hadError.set(true);
            } else {
                PointValuePair result = resultOpt.get();
                double[] rPoint = result.getPoint();
                for (int j = 0; j < nPar; j++) {
                    parValues[j][i] = rPoint[j];
                }
                parValues[nPar][i] = result.getValue();

                if (equation == CPMGEquation.CPMGSLOW) {
                    for (int j = 0; j < map.length; j++) {
                        rexValues[j][i] = equation.getRex(result.getPoint(), map[j], xValues[1][0]);
                    }
                }
            }
        });

        if (hadError.get()) {
            return Optional.empty();
        } else {
            double[] parSDev = new double[nPar];
            for (int i = 0; i < nPar; i++) {
                DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
                double v25 = dStat.getPercentile(25);
                double v75 = dStat.getPercentile(75.0);
                parSDev[i] = (v75 - v25) / 1.349;
            }
            if (equation == CPMGEquation.CPMGSLOW) {
                for (int j = 0; j < nID; j++) {
                    DescriptiveStatistics dStat = new DescriptiveStatistics(rexValues[j]);
                    rexErrors[j] = dStat.getStandardDeviation();
                }
            }
            return Optional.of(parSDev);
        }
    }

}
