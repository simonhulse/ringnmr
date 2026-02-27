package org.comdnmr.eqnfit;

import org.apache.commons.math3.optim.PointValuePair;
import org.checkerframework.checker.nullness.Opt;
import org.comdnmr.fit.FitQuality;
import org.comdnmr.util.CoMDOptions;
import org.nmrfx.chemistry.relax.ResonanceSource;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class SSR1RhoFitter implements EquationFitter {

    FitFunction fitFunc;
    CoMDOptions options;

    List<Double>[] xValues;
    List<Double> yValues = new ArrayList<>();
    List<Double> errValues = new ArrayList<>();
    List<Integer> idValues = new ArrayList<>();
    int nCurves = 1;
    int nResidues = 1;
    int[][] states;
    int[] stateCount;
    ResonanceSource[] dynSources;
    static List<String> equationNameList = Arrays.asList(SSR1RhoEquation.getAllEquationNames());
    long errTime;
    static final String expType = "ssr1rho";

    public SSR1RhoFitter(CoMDOptions options) {
        fitFunc = new SSR1RhoFitFunction(options);
        this.options = options;
    }

    @Override
    public List<String> getEquationNameList() {
        return Arrays.asList(SSR1RhoEquation.getAllEquationNames());
    }

    @Override
    public FitFunction getFitModel() {
        return fitFunc;
    }

    @Override
    public String getExpType() {
        return "";
    }

    @Override
    public Optional<FitResult> doFit(String eqn, double[] sliderGuesses, CoMDOptions options) {
        setupFit(eqn);
        int[][] map = fitFunc.getMap();
        double[] guesses;
        if (sliderGuesses != null) {
            //fixme
            guesses = sliderGuesses;
        } else {
            guesses = fitFunc.guess();
        }
        double[][] boundaries = fitFunc.boundaries(guesses);
        double sigma = options.getStartRadius();
        var resultOpt = fitFunc.refine(guesses, boundaries[0], boundaries[1], sigma, options.getOptimizer());
        if (resultOpt.isEmpty()) {
            return Optional.empty();
        }
        PointValuePair result = resultOpt.get();
        double[] pars = result.getPoint();
        FitQuality fitQuality = fitFunc.getFitQuality(pars);
        int nGroupPars = fitFunc.getNGroupPars();
        sigma /= 2.0;

        String[] parNames = fitFunc.getParNames();
        double[] errEstimates = new double[pars.length];
        double[][] simPars = null;
        if (FitFunction.getCalcError()) {
            long startTime = System.currentTimeMillis();
            var errOpt = fitFunc.simBoundsStream(pars.clone(),
                    boundaries[0], boundaries[1], sigma, options);
            long endTime = System.currentTimeMillis();
            if (errOpt.isPresent()) {
                errEstimates = errOpt.get();
            }

            errTime = endTime - startTime;
            simPars = fitFunc.getSimPars();
        }
        String refineOpt = options.getOptimizer();
        String bootstrapOpt = options.getBootStrapOptimizer();
        long fitTime = fitFunc.fitTime;
        long bootTime = errTime;
        int nSamples = options.getSampleSize();
        boolean useAbs = options.getAbsValueFit();
        boolean useNonParametric = options.getNonParametricBootstrap();
        double sRadius = options.getStartRadius();
        double fRadius = options.getFinalRadius();
        double tol = options.getTolerance();
        boolean useWeight = options.getWeightFit();
        CurveFit.CurveFitStats curveStats = new CurveFit.CurveFitStats(refineOpt, bootstrapOpt, fitTime, bootTime, nSamples, useAbs,
                useNonParametric, sRadius, fRadius, tol, useWeight);
        double[][] extras = getFields(xValues, idValues);
        FitResult res = getResults(this, eqn, parNames, dynSources, map, states, extras, nGroupPars, pars, errEstimates, fitQuality, simPars, true, curveStats);
        return Optional.of(res);
    }

    @Override
    public int[] getStateCount() {
        return new int[0];
    }

    @Override
    public int[][] getStates() {
        return new int[0][];
    }

    // NOTE: Copied over from R1RhoFitter
    // TODO: Could make a utility method called linespace to do this in one line.
    // Probably useful elsewhere in the codebase too
    @Override
    public double[] getSimX(int nPts, double xLB, double xUB) {
        int nPoints = nPts;
        double[] x = new double[nPoints];
        double firstValue = xLB;
        double lastValue = xUB;
        double delta = (lastValue - firstValue) / (nPoints + 1);
        double value = firstValue;
        for (int i = 0; i < nPoints; i++) {
            x[i] = value;
            value += delta;
        }
        return x;
    }

    @Override
    public double rms(double[] pars) {
        return 0;
    }

    @Override
    public void setData(List<Double>[] allXValues, List<Double> yValuesLog10, List<Double> errValuesLog10) {
        xValues = new ArrayList[allXValues.length];
        for (int j = 0; j < allXValues.length; j++) {
            xValues[j] = new ArrayList<>();
            xValues[j].addAll(allXValues[j]);
        }
        this.yValues.clear();

        // Y-axis for SS R1rho is logarithmic.
        // Need to convert log10(y) -> y
        // and log10(err_y) -> err_y
        List<Double> yValues = new ArrayList<>(yValuesLog10.size());
        List<Double> errValues = new ArrayList<>(errValuesLog10.size());
        for (int i = 0; i < yValuesLog10.size(); i++) {
            yValues.add(Math.pow(10.0, yValuesLog10.get(i)));
            errValues.add(Math.pow(10.0, errValuesLog10.get(i)));
        }
        this.yValues.addAll(yValues);
        this.errValues.clear();
        this.errValues.addAll(errValues);
        this.idValues.clear();
        yValues.forEach((_item) -> {
            this.idValues.add(0);
        });
        dynSources = new ResonanceSource[1];
        dynSources[0] = null;
        nCurves = 1;
        stateCount = new int[4];
        stateCount[0] = nResidues;
        stateCount[1] = 1;
        stateCount[2] = 1;
        stateCount[3] = 1;
        states = new int[1][4];
    }

    @Override
    public void setupFit(String eqn) {
        double[][] x = new double[xValues.length][yValues.size()];
        double[] y = new double[yValues.size()];
        double[] err = new double[yValues.size()];
        int[] idNums = new int[yValues.size()];
        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0; j < x.length; j++) {
                x[j][i] = xValues[j].get(i);
            }
            y[i] = yValues.get(i);
            err[i] = errValues.get(i);
            idNums[i] = idValues.get(i);
        }
        fitFunc.setEquation(eqn);
        fitFunc.setXY(x, y);
        fitFunc.setIds(idNums);
        fitFunc.setErr(err);
        fitFunc.setMap(stateCount, states);
    }


    @Override
    public List<ParValueInterface> guessPars(String eqn) {
        return List.of();
    }

    @Override
    public double[] getSimXDefaults() {
        return new double[0];
    }
}
