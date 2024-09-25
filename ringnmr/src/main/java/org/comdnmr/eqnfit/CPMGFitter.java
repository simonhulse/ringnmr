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

import org.comdnmr.data.CPMGExperiment;
import org.comdnmr.data.ExperimentSet;
import org.comdnmr.data.ExperimentData;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.inference.TTest;
import org.comdnmr.data.Experiment;
import org.comdnmr.fit.FitQuality;
import org.comdnmr.util.CoMDOptions;
import org.comdnmr.util.CoMDPreferences;
import org.nmrfx.chemistry.relax.ResonanceSource;

/**
 *
 * @author Bruce Johnson
 */
public class CPMGFitter implements EquationFitter {

    public static final double[] SIMX = {25.0, 50.0, 100.0, 150.0, 250.0, 350.0, 500.0, 750.0, 1000.0};

    CPMGFitFunction calcR;
    CoMDOptions options;
    List<Double>[] xValues = new ArrayList[4];
    List<Double> yValues = new ArrayList<>();
    List<Double> errValues = new ArrayList<>();
    List<Integer> idValues = new ArrayList<>();
    int nCurves = 1;
    int nResidues = 1;
    int[][] states;
    int[] stateCount;
    ResonanceSource[] dynSources;
    long errTime;
    static final String expType = "cpmg";

    class StateCount {

        int[][] states;
        int[] stateCount;

        StateCount(int[][] states, int[] stateCount) {
            this.states = states;
            this.stateCount = stateCount;
        }

        int getResIndex(int i) {
            return states[i][0];
        }

        int getTempIndex(int i) {
            return states[i][2];
        }
    }

    public CPMGFitter(CoMDOptions options) {
        calcR = new CPMGFitFunction(options);
        this.options = options;
    }

    public static int getMapIndex(int[] state, int[] stateCount, int... mask) {
        int index = 0;
        double mult = 1.0;
        for (int j : mask) {
            index += mult * state[j];
            mult *= stateCount[j];
        }
        return index;
    }

    @Override
    public void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues) {
        for (int i = 0;i<allXValues.length;i++) {
            xValues[i] = new ArrayList<>();
            xValues[i].addAll(allXValues[i]);
        }
        this.yValues.clear();
        this.yValues.addAll(yValues);
        this.errValues.clear();
        this.errValues.addAll(errValues);
        this.idValues.clear();
        for (Double yValue : yValues) {
            this.idValues.add(0);
        }
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
    public void setData(ExperimentSet experimentSet, ResonanceSource[] dynSources) {
        for (int i = 0;i<xValues.length;i++) {
            xValues[i] = new ArrayList<>();
        }
        this.dynSources = dynSources.clone();
        nResidues = dynSources.length;

        stateCount = experimentSet.getStateCount(nResidues);
        Collection<Experiment> expDataList = experimentSet.getExperimentData();
        nCurves = dynSources.length * expDataList.size();
        states = new int[nCurves][];
        int k = 0;
        int resIndex = 0;
        int id = 0;
        for (var dynSource : dynSources) {
            for (Experiment experiment : expDataList) {
                CPMGExperiment cpmgExperiment = (CPMGExperiment) experiment;
                ExperimentData experimentalData = experiment.getResidueData(dynSource);
                if (experimentalData != null) {

                    states[k++] = experimentSet.getStateIndices(resIndex, experiment);
                    //  need peakRefs
                    double fieldX = experiment.getNucleusField();
                    double fieldH = experiment.getB0Field();
                    double[][] x = experimentalData.getXValues();
                    double[] y = experimentalData.getYValues();
                    double[] err = experimentalData.getErrValues();

                    for (int i = 0; i < y.length; i++) {
                        xValues[0].add(x[0][i]);
                        xValues[1].add(fieldX);
                        xValues[2].add(fieldH);
                        xValues[3].add(cpmgExperiment.getTau());
                        yValues.add(y[i]);
                        errValues.add(err[i]);
                        idValues.add(id);
                    }
                    id++;
                }
            }
            resIndex++;
        }
    }

    @Override
    public FitFunction getFitModel() {
        return calcR;
    }

    @Override
    public List<String> getEquationNameList() {
        return getEquationNames();
    }
    
    @Override
    public String getExpType() {
        return expType;
    }

    public static List<String> getEquationNames() {
        List<String> activeEquations = CoMDPreferences.getActiveCPMGEquations();
        System.out.println(activeEquations);
        return activeEquations;
    }

    @Override
    public int[] getStateCount() {
        return stateCount;
    }

    @Override
    public int[][] getStates() {
        return states;
    }

    @Override
    public void setupFit(String eqn) {
        double[][] x = new double[xValues.length][yValues.size()];
        double[] y = new double[yValues.size()];
        double[] err = new double[yValues.size()];
        int[] idNums = new int[yValues.size()];
        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0;j < xValues.length;j++) {
                if (xValues[j] != null) {
                    x[j][i] = xValues[j].get(i);
                }
            }
            y[i] = yValues.get(i);
            err[i] = errValues.get(i);
            idNums[i] = idValues.get(i);
        }
        calcR.setEquation(eqn);

        calcR.setXY(x, y);
        calcR.setIds(idNums);
        calcR.setErr(err);
        calcR.setMap(stateCount, states);
    }

    @Override
    public List<ParValueInterface> guessPars(String eqn) {
        setupFit(eqn);
        double[] guesses = calcR.guess();
        String[] parNames = calcR.getParNames();
        int[][] map = calcR.getMap();
        List<ParValueInterface> parValues = new ArrayList<>();
        for (int i = 0; i < parNames.length; i++) {
            double guess = guesses[map[0][i]];
            ParValueInterface parValue = new ParValue(parNames[i], guess);
            parValues.add(parValue);
        }
        return parValues;
    }

    @Override
    public double rms(double[] pars) {
        return calcR.getRMS(pars);
    }

    @Override
    public FitResult doFit(String eqn, double[] sliderguesses, CoMDOptions options) {
        setupFit(eqn);
        int[][] map = calcR.getMap();
        double[] guesses;
        if (sliderguesses != null) {
            //fixme
            guesses = sliderguesses;
        } else {
            guesses = calcR.guess();
        }
        double[][] boundaries = calcR.boundaries(guesses);
        double sigma = options.getStartRadius();
        PointValuePair result = calcR.refine(guesses, boundaries[0],
                boundaries[1], sigma, options.getOptimizer());
        double[] pars = result.getPoint();
        System.out.print("Fit pars \n");
        for (int i = 0; i < pars.length; i++) {
            System.out.printf("%d %.3f %.3f %.3f %.3f\n", i, guesses[i], boundaries[0][i], pars[i], boundaries[1][i]);
        }
        System.out.println();

        FitQuality fitQuality = calcR.getFitQuality(pars);

        int nGroupPars = calcR.getNGroupPars();
        sigma /= 2.0;

        String[] parNames = calcR.getParNames();
        double[] errEstimates;
        double[][] simPars = null;
        double field = calcR.xValues[1][0];
        double[] rexValues = calcR.getRex(pars, field);
        boolean okRex = false;
        double rexRatio = options.getRexRatio();
        for (double rexValue : rexValues) {
            if (rexValue > rexRatio * fitQuality.rms()) {
                okRex = true;
                break;
            }
        }
        boolean exchangeValid = okRex;

        if (FitFunction.getCalcError()) {
            long startTime = System.currentTimeMillis();
            errEstimates = calcR.simBoundsStream(pars.clone(),
                    boundaries[0], boundaries[1], sigma, options);
            long endTime = System.currentTimeMillis();
            errTime = endTime - startTime;
            simPars = calcR.getSimPars();
            TTest tTest = new TTest();
            for (String parName : parNames) {
                int parIndex = -1;
                if (parName.equals("Kex")) {
                    parIndex = 0;
                    if (pars[parIndex] < errEstimates[parIndex]) {
                        exchangeValid = false;
                    }
                }
                if (parName.equals("dPPMmin")) {
                    parIndex = 2;
                }
                if (parIndex != -1) {
                    boolean valid = tTest.tTest(0.0, simPars[map[0][parIndex]], 0.02);
                    SummaryStatistics sStat = new SummaryStatistics();
                    for (double v : simPars[map[0][parIndex]]) {
                        sStat.addValue(v);
                    }
                    double alpha = tTest.tTest(0.0, simPars[map[0][parIndex]]);
                    double mean = sStat.getMean();
                    double sdev = sStat.getStandardDeviation();
                    if (!valid) {
                        exchangeValid = false;
                    }
                }
            }
        } else {
            errEstimates = new double[pars.length];
        }
        String refineOpt = options.getOptimizer();
        String bootstrapOpt = options.getBootStrapOptimizer();
        long fitTime = calcR.fitTime;
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
        return getResults(this, eqn, parNames, dynSources, map, states, extras, nGroupPars, pars, errEstimates, fitQuality, simPars, exchangeValid, curveStats);
    }

    @Override
    public double[] getSimX(int nPts, double xLB, double xUB) {
        int nPoints = nPts;
        double[] x = new double[nPoints];
        double delta = (xUB - xLB) / (nPoints + 1);
        double value = xLB;
        for (int i = 0; i < nPoints; i++) {
            x[i] = value;
            value += delta;

        }
        return x;
    }

    @Override
    public double[] getSimXDefaults() {
        return SIMX;
    }
}
