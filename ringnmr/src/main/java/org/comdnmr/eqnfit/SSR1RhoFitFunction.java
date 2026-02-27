package org.comdnmr.eqnfit;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;
import org.comdnmr.util.CoMDOptions;

import java.util.Optional;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;

public class SSR1RhoFitFunction extends FitFunction {

    static RandomGenerator random = new SynchronizedRandomGenerator(new Well19937c());

    public SSR1RhoFitFunction(CoMDOptions options) {
        super(options);
        this.equation = SSR1RhoEquation.CSA;
    }

    public SSR1RhoFitFunction(CoMDOptions options, double[][] x, double[] y, double[] err) throws IllegalArgumentException {
        this(options, x, y, err, new int[x.length]);
    }

    public SSR1RhoFitFunction(CoMDOptions options, double[][] x, double[] y, double[] err, int[] idNums) throws IllegalArgumentException {
        super(options);
        this.xValues = new double[1][];
        this.xValues[0] = x[0].clone();
        this.yValues = y.clone();
        this.errValues = err.clone();
        this.idNums = idNums.clone();
        this.idNums = new int[yValues.length];
        this.equation = SSR1RhoEquation.CSA;
        if (setNID()) {
            throw new IllegalArgumentException("Invalid idNums, some values not used");
        }
    }

    @Override
    public void setEquation(String eqName) {
        equation = SSR1RhoEquation.valueOf(eqName.toUpperCase());
    }

    // TODO
    @Override
    public int[] getMask() {
        return new int[0];
    }

    // TODO: simBounds has no usages. Delete from FitFunction
    @Override
    public Optional<double[]> simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        return Optional.empty();
    }

    @Override
    public double[][] getSimPars() { return parValues; }

    // Currently only set up for simple case, without any mappings/IDs
    @Override
    public double value(double[] normPar) {
        double[] par = deNormalize(normPar);
        double[][] x = xValues;
        double[] ys = equation.calculate(par, map[0], x, 0);

        double result = 0.0;
        for (int i = 0; i < ys.length; i++) {
            double delta = ys[i] - yValues[i];
            if (weightFit) delta /= errValues[i];
            if (absMode) {
                result += FastMath.abs(delta);
            } else {
                result += FastMath.pow(delta, 2);
            }
        }
        return result / (yValues.length - par.length);
    }

    // NOTE:
    // simBoundsStream, simBoundsStreamNonParametric, and simBoundsStreamParametric
    // are all copy-pasted from ExpFitFunction
    @Override
    public Optional<double[]> simBoundsStream(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        if (options.getNonParametricBootstrap()) {
            return simBoundsStreamNonParametric(start, lowerBounds, upperBounds, inputSigma, options);
        } else {
            return simBoundsStreamParametric(start, lowerBounds, upperBounds, inputSigma, options);
        }
    }

    public Optional<double[]> simBoundsStreamParametric(double[] start,
                                                        double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = options.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        double[] yPred = getPredicted(start);
        String optimizer = options.getBootStrapOptimizer();
        AtomicBoolean failed = new AtomicBoolean(false);

        IntStream.range(0, nSim).parallel().forEach(i -> {
//        IntStream.range(0, nSim).forEach(i -> {
            ExpFitFunction rDisp = new ExpFitFunction(options, xValues, yPred, errValues, idNums);
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
                failed.set(true);
            } else {
                PointValuePair result = resultOpt.get();
                double[] rPoint = result.getPoint();
                for (int j = 0; j < nPar; j++) {
                    parValues[j][i] = rPoint[j];
                }
                parValues[nPar][i] = result.getValue();
            }
        });

        if (failed.get()) {
            return Optional.empty();
        }
        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        return Optional.of(parSDev);
    }

    public Optional<double[]> simBoundsStreamNonParametric(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = options.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        String optimizer = options.getBootStrapOptimizer();
        AtomicBoolean failed = new AtomicBoolean(false);

        IntStream.range(0, nSim).parallel().forEach(i -> {
            SSR1RhoFitFunction rDisp = new SSR1RhoFitFunction(options, xValues, yValues, errValues, idNums);
            rDisp.setEquation(equation.getName());
            double[][] newX = new double[xValues.length][yValues.length];
            double[] newY = new double[yValues.length];
            double[] newErr = new double[yValues.length];
            int[] newID = new int[yValues.length];
            int iTry = 0;
            do {
                for (int j = 0; j < xValues.length; j++) {
                    for (int k = 0; k < yValues.length; k++) {
                        int rI = random.nextInt(yValues.length);
                        newX[j][k] = xValues[j][rI];
                        newY[k] = yValues[rI];
                        newErr[k] = errValues[rI];
                        newID[k] = idNums[rI];
                    }
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
                failed.set(true);
            } else {
                PointValuePair result = resultOpt.get();
                double[] rPoint = result.getPoint();
                for (int j = 0; j < nPar; j++) {
                    parValues[j][i] = rPoint[j];
                }
                parValues[nPar][i] = result.getValue();
            }
        });
        if (failed.get()) {
            return Optional.empty();
        }

        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        return Optional.of(parSDev);
    }
}
