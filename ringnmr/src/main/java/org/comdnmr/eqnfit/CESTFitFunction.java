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

public class CESTFitFunction extends FitFunction {

    static RandomGenerator random = new SynchronizedRandomGenerator(new Well19937c());

    int[] r2Mask = {0, 1, 3};
    double[] rexErrors = new double[nID];

    public CESTFitFunction(CoMDOptions options) {
        super(options);
        this.equation = CESTEquation.TROTT_PALMER;
    }

    @Override
    public void setEquation(String eqName) {
        equation = CESTEquation.valueOf(eqName.toUpperCase());
    }

    public CESTFitFunction(CoMDOptions options, double[][] x, double[] y, double[] err) throws IllegalArgumentException {
        this(options, x, y, err, new int[x.length]);
    }

    public CESTFitFunction(CoMDOptions options, double[][] x, double[] y, double[] err, int[] idNums) throws IllegalArgumentException {
        super(options);
        this.xValues = new double[x.length][];
        for (int j = 0; j < x.length; j++) {
            this.xValues[j] = x[j].clone();
        }
        this.yValues = y.clone();
        this.errValues = err.clone();
        this.idNums = idNums.clone();
        this.idNums = new int[yValues.length];
        this.equation = CESTEquation.TROTT_PALMER;
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
    public double value(double[] normPar) {
        double[] par = deNormalize(normPar);

        double sumAbs = 0.0;
        double sumSq = 0.0;
        double[] yCalc = new double[yValues.length];
        for (int id = 0; id < map.length; id++) {
            double[][] x = CESTEquations.getXValues(xValues, idNums, id);
            double[] yCalc1 = equation.calculate(par, map[id], x, id);
            int[] indicies = CESTEquations.getIndicies(idNums, id);
            for (int i = 0; i < indicies.length; i++) {
                yCalc[indicies[i]] = yCalc1[i];
            }
        }

        for (int i = 0; i < yValues.length; i++) {
            double delta = (yCalc[i] - yValues[i]);
            if (weightFit) {
                delta /= errValues[i];
            }
            sumAbs += FastMath.abs(delta);
            sumSq += delta * delta;
        }
        if (absMode) {
            return sumAbs / (yValues.length - par.length);
        } else {
            return sumSq / (yValues.length - par.length);
        }
    }

    @Override
    public double[] getPredicted(double[] par) {
        return simY(par);
    }

    @Override
    public double[][] getSimPars() {
        return parValues;
    }

    @Override
    public Optional<double[]> simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = options.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        double[] yPred = getPredicted(start);
        double[] yValuesOrig = yValues.clone();
        rexErrors = new double[nID];
        String optimizer = options.getBootStrapOptimizer();
        for (int i = 0; i < nSim; i++) {
            for (int k = 0; k < yValues.length; k++) {
                yValues[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            var resultOpt = refine(start, lowerBounds, upperBounds,
                    inputSigma, optimizer);
            if (resultOpt.isEmpty()) {
                return Optional.empty();
            }
            PointValuePair result = resultOpt.get();
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
            parValues[nPar][i] = result.getValue();

        }
        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        yValues = yValuesOrig;
        return Optional.of(parSDev);
    }

    @Override
    public Optional<double[]> simBoundsStream(double[] start, double[] lowerBounds,
                                              double[] upperBounds, double inputSigma, CoMDOptions options) {
        if (Boolean.TRUE.equals(options.getNonParametricBootstrap())) {
            return simBoundsStreamNonParametric(start, lowerBounds, upperBounds, inputSigma, options);
        } else {
            return simBoundsStreamParametric(start, lowerBounds, upperBounds, inputSigma, options);
        }

    }

    public Optional<double[]> simBoundsStreamParametric(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = options.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        rexErrors = new double[nID];
        double[] yPred = getPredicted(start);
        String optimizer = options.getBootStrapOptimizer();
        AtomicBoolean hadError = new AtomicBoolean(false);

        IntStream.range(0, nSim).parallel().forEach(i -> {
            CESTFitFunction rDisp = new CESTFitFunction(options, xValues, yPred, errValues, idNums);
            rDisp.setEquation(equation.getName());
            double[] newY = new double[yValues.length];
            for (int k = 0; k < yValues.length; k++) {
                newY[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            rDisp.setXY(xValues, newY);
            rDisp.setIds(idNums);
            rDisp.setMap(map);

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
            }
        });
        if (hadError.get()) {
            return Optional.empty();
        } else {
            double[] parSDev = new double[nPar];
            for (int i = 0; i < nPar; i++) {
                DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
                parSDev[i] = dStat.getStandardDeviation();
            }
            return Optional.of(parSDev);
        }
    }

    public Optional<double[]> simBoundsStreamNonParametric(double[] start,
                                                           double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = options.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        rexErrors = new double[nID];
        String optimizer = options.getBootStrapOptimizer();

        AtomicBoolean hadError = new AtomicBoolean(false);
        IntStream.range(0, nSim).parallel().forEach(i -> {
            CESTFitFunction rDisp = new CESTFitFunction(options, xValues, yValues, errValues, idNums);
            rDisp.setEquation(equation.getName());
            double[][] newX = new double[xValues.length][yValues.length];
            double[] newY = new double[yValues.length];
            double[] newErr = new double[yValues.length];
            int[] newID = new int[yValues.length];
            int iTry = 0;
            do {
                for (int k = 0; k < yValues.length; k++) {
                    int rI = random.nextInt(yValues.length);
                    for (int j = 0; j < xValues.length; j++) {
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
            }
        });
        if (hadError.get()) {
            return Optional.empty();
        } else {
            double[] parSDev = new double[nPar];
            for (int i = 0; i < nPar; i++) {
                DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
                parSDev[i] = dStat.getStandardDeviation();
            }
            return Optional.of(parSDev);
        }
    }

}
