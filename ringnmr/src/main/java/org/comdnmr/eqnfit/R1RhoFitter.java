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

import org.comdnmr.fit.FitQuality;
import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.data.ExperimentSet;
import org.comdnmr.data.ExperimentData;
import org.comdnmr.data.Experiment;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.optim.PointValuePair;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.comdnmr.util.CoMDOptions;
import org.nmrfx.chemistry.Atom;

/**
 *
 * @author Bruce Johnson
 */
public class R1RhoFitter implements EquationFitter {

    FitFunction calcR1Rho;
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
    static List<String> equationNameList = Arrays.asList(R1RhoEquation.getEquationNames());
    long errTime;
    Map<String, List<Double>>[] constraints;
    static final String expType = "r1rho";

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

    public R1RhoFitter(CoMDOptions options) {
        calcR1Rho = new R1RhoFitFunction(options);
        this.options = options;
    }

    public static int getMapIndex(int[] state, int[] stateCount, int... mask) {
        int index = 0;
//        System.out.println(state.length + " mask " + mask.length);
//        for (int i = 0; i < state.length; i++) {
//            System.out.print(" " + state[i]);
//        }
//        System.out.println("");
        double mult = 1.0;
        for (int i = 0; i < mask.length; i++) {
//            System.out.println("mask:" + mask[i] + " state[mask]:" + state[mask[i]] + " count:" + stateCount[mask[i]]);
            index += mult * state[mask[i]];
            mult *= stateCount[mask[i]];
        }
        return index;
    }

    @Override
    public void setData(ExperimentSet experimentSet, ResonanceSource[] dynSources) {
        xValues = new ArrayList[4];
        for (int j = 0;j<xValues.length;j++) {
            xValues[j] = new ArrayList<>();
        }
        this.dynSources = dynSources.clone();
        nResidues = dynSources.length;
        experimentSet.setupMaps();
        stateCount = experimentSet.getStateCount(dynSources.length);
        Collection<Experiment> expDataList = experimentSet.getExperimentData();
        nCurves = experimentSet.getDataCount(dynSources);
        states = new int[nCurves][];
        int k = 0;
        int resIndex = 0;
        int id = 0;
        constraints = new Map[nCurves];
        for (var atom : dynSources) {
            for (Experiment expData : expDataList) {
                ExperimentData experimentalData = expData.getResidueData(atom);
                if (experimentalData != null) {
                    constraints[id] = expData.getConstraints();
                    states[k++] = experimentSet.getStateIndices(resIndex, expData);
                    double[][] x = experimentalData.getXValues();
                    double[] y = experimentalData.getYValues();
                    double[] err = experimentalData.getErrValues();
                    double field = expData.getNucleusField();
                    for (int i = 0; i < y.length; i++) {
                        for (int j=0;j<x.length;j++) {
                            xValues[j].add(x[j][i]);
                        }
                        xValues[xValues.length - 1].add(field);
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
    public void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues) {
        xValues = new ArrayList[allXValues.length];
        for (int j=0;j<allXValues.length;j++) {
            xValues[j] = new ArrayList<>();
            xValues[j].addAll(allXValues[j]);
        }
        this.yValues.addAll(yValues);
        this.errValues.addAll(errValues);
        for (Double yValue : yValues) {
            idValues.add(0);
        }
        dynSources = new ResonanceSource[1];
        dynSources[0] = null;

        stateCount = new int[4];
        stateCount[0] = 1;
        stateCount[1] = 1;
        stateCount[2] = 1;
        stateCount[3] = 1;

        states = new int[1][4];
        states[0][0] = 0;
        states[0][1] = 0;
        states[0][2] = 0;
        states[0][3] = 0;
    }

    @Override
    public FitFunction getFitModel() {
        return calcR1Rho;
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
        List<String> activeEquations = CoMDPreferences.getActiveR1RhoEquations();
        System.out.println(activeEquations.toString());
        return activeEquations;

    }

    public static List<String> setEquationNames(List<String> eqnNames) {
        equationNameList = eqnNames;
        return equationNameList;
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
            for (int j=0;j<xValues.length;j++) {
                x[j][i] = xValues[j].get(i);
            }
            y[i] = yValues.get(i);
            err[i] = errValues.get(i);
            idNums[i] = idValues.get(i);
        }
        calcR1Rho.setEquation(eqn);
        calcR1Rho.setXY(x, y);
        calcR1Rho.setIds(idNums);
        calcR1Rho.setErr(err);
        calcR1Rho.setMap(stateCount, states);
    }

    @Override
    public List<ParValueInterface> guessPars(String eqn) {
        setupFit(eqn);
        double[] guesses = calcR1Rho.guess();
        if (guesses != null) {
            String[] parNames = calcR1Rho.getParNames();
            List<ParValueInterface> parValues = new ArrayList<>();
            int[][] map = calcR1Rho.getMap();
            for (int i = 0; i < parNames.length; i++) {
                ParValueInterface parValue = new ParValue(parNames[i], guesses[map[0][i]]);
                parValues.add(parValue);
            }
            return parValues;
        } else {
            return null;
        }
    }

    @Override
    public double rms(double[] pars) {
        double rms = calcR1Rho.getRMS(pars);
        return rms;
    }

    @Override
    public FitResult doFit(String eqn, double[] sliderguesses, CoMDOptions options) {
        double[][] xvals = new double[xValues.length][xValues[0].size()];
        double[] yvals = new double[yValues.size()];
        int[] idNums = new int[yValues.size()];
        for (int i = 0; i < xvals.length; i++) {
            for (int j = 0; j < xvals[0].length; j++) {
                xvals[i][j] = xValues[i].get(j);
            }
        }
        for (int i = 0; i < yvals.length; i++) {
            yvals[i] = yValues.get(i);
            idNums[i] = idValues.get(i);
        }
        double[][] xy = CESTEquations.getXYValues(xvals, yvals, idNums, 0);
        List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy, "r1rho");
        if (peaks.size() >= 1) {
            setupFit(eqn);
            int[][] map = calcR1Rho.getMap();
            double[] guesses;
            if (sliderguesses != null) {
                //fixme
                guesses = sliderguesses;
            } else {
                guesses = calcR1Rho.guess();
            }
            if (guesses != null) {
                double[][] boundaries = calcR1Rho.boundaries(guesses);
                double sigma = options.getStartRadius();
                if (constraints != null) {
                    for (int id = 0; id < map.length; id++) {
                        if (constraints[id] != null) {
                            List<Double> cValues = constraints[id].get("R1A");
                            if (cValues != null) {
                                double[] cArray = {cValues.get(0), cValues.get(1)};
                                calcR1Rho.equation.constrain("R1A", guesses, boundaries, map, id, cArray[0], cArray[1]);
                                calcR1Rho.equation.constrain("R1B", guesses, boundaries, map, id, cArray[0], cArray[1]);
                            }
                        }
                    }
                }

                for (int i = 0; i < guesses.length; i++) {
                    System.out.println(i + " bou0 " + boundaries[0][i] + " bou1 " + boundaries[1][i] + " gue " + guesses[i]);
                }
                PointValuePair result = calcR1Rho.refine(guesses, boundaries[0],
                        boundaries[1], sigma, options.getOptimizer());
                double[] pars = result.getPoint();
                System.out.println(eqn);

                for (int[] map1 : map) {
                    for (int j = 0; j < map1.length; j++) {
                        System.out.printf(" %3d", map1[j]);
                    }
                    System.out.println("");
                }

                System.out.print("Fit pars ");
                for (int i = 0; i < pars.length; i++) {
                    System.out.printf(" %.3f", pars[i]);
                }
                System.out.println("");

                FitQuality fitQuality = calcR1Rho.getFitQuality(pars);

                int nGroupPars = calcR1Rho.getNGroupPars();
                sigma /= 2.0;

                String[] parNames = calcR1Rho.getParNames();
                double[] errEstimates;
                double[][] simPars = null;
                boolean exchangeValid = true;
                double deltaABdiff = options.getDeltaABDiff();
                if (FitFunction.getCalcError()) {
                    long startTime = System.currentTimeMillis();
                    errEstimates = calcR1Rho.simBoundsStream(pars.clone(),
                            boundaries[0], boundaries[1], sigma, options);
                    long endTime = System.currentTimeMillis();
                    errTime = endTime - startTime;

                    simPars = calcR1Rho.getSimPars();
                    for (String parName : parNames) {
                        if (parName.equals("deltaB0")) {
                            int parIndex = 3;
                            int deltaAIndex = parIndex - 1;
                            if (Math.abs(pars[parIndex] - pars[deltaAIndex]) < deltaABdiff) {
                                exchangeValid = false;
                            }
                        }
                    }
                } else {
                    errEstimates = new double[pars.length];
                }
                double[][] extras = getFields(xValues, idValues);
                String refineOpt = options.getOptimizer();
                String bootstrapOpt = options.getBootStrapOptimizer();
                long fitTime = calcR1Rho.fitTime;
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
                return getResults(this, eqn, parNames, dynSources, map, states, extras, nGroupPars, pars, errEstimates, fitQuality, simPars, exchangeValid, curveStats);
            } else {
                return null;
            }
        } else {
            return null;
        }
    }

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
    public double[] getSimXDefaults() {
        return getSimX(100, -8, 8);
    }

}
