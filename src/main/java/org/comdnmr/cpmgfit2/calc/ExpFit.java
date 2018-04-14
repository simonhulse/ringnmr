package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.math3.optim.PointValuePair;

/**
 *
 * @author Bruce Johnson
 */
public class ExpFit implements EquationFitter {

    List<Double> xValues = new ArrayList<>();
    List<Double> yValues = new ArrayList<>();
    List<Double> errValues = new ArrayList<>();
    List<Double> fieldValues = new ArrayList<>();
    List<Integer> idValues = new ArrayList<>();
    double[] usedFields = null;
    int nCurves = 1;
    int nResidues = 1;
    int[][] states;
    int[] stateCount;
    String[] resNums;
    static List<String> equationNameList = Arrays.asList(ExpEquation.getEquationNames());

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

    // public void setData(Collection<ExperimentData> expDataList, String[] resNums) {
    public void setData(ResidueProperties resProps, String[] resNums) {
        this.resNums = resNums.clone();
        nResidues = resNums.length;
        int id = 0;
        resProps.setupMaps();
        stateCount = resProps.getStateCount(resNums.length);
        Collection<ExperimentData> expDataList = resProps.getExperimentData();
        nCurves = resNums.length * expDataList.size();
        states = new int[nCurves][];
        int k = 0;
        int resIndex = 0;
        for (String resNum : resNums) {
            for (ExperimentData expData : expDataList) {
                states[k++] = resProps.getStateIndices(resIndex, expData);
                ResidueData resData = expData.getResidueData(resNum);
                //  need peakRefs
                double field = expData.getField();
                double[] x = resData.getXValues();
                double[] y = resData.getYValues();
                double[] err = resData.getErrValues();
                for (int i = 0; i < x.length; i++) {
                    xValues.add(x[i]);
                    yValues.add(y[i]);
                    errValues.add(err[i]);
                    fieldValues.add(field);
                    idValues.add(id);
                }
                id++;

            }
            resIndex++;
        }
        usedFields = new double[expDataList.size()];
        int iExp = 0;
        for (ExperimentData expData : expDataList) {
            usedFields[iExp++] = expData.getField();
        }
    }

    public static List<String> getEquationNames() {
        return equationNameList;
    }

    public CPMGFitResult doFit(ResidueProperties resProps, String eqn, boolean absMode, boolean nonParBootStrap) {
        double[] x = new double[xValues.size()];
        double[] y = new double[xValues.size()];
        double[] err = new double[xValues.size()];
        int[] idNums = new int[xValues.size()];
        double[] fields = new double[xValues.size()];
        for (int i = 0; i < x.length; i++) {
            x[i] = xValues.get(i);
            y[i] = yValues.get(i);
            err[i] = errValues.get(i);
            fields[i] = fieldValues.get(i);
            idNums[i] = idValues.get(i);
        }
        FitModel calcR = new CalcExpDecay();
        calcR.setEquation(eqn);
        calcR.setAbsMode(absMode);

        calcR.setXY(x, y);
        calcR.setIds(idNums);
        calcR.setErr(err);
        calcR.setFieldValues(fields);
        calcR.setFields(usedFields);
        calcR.setMap(stateCount, states);
        int[][] map = calcR.getMap();
        double[] guesses = calcR.guess();
        double[][] boundaries = calcR.boundaries();
        double[] sigma = new double[guesses.length];
        for (int i = 0; i < guesses.length; i++) {
            sigma[i] = (boundaries[1][i] - boundaries[0][i]) / 10.0;
//            System.out.println(i + " " + boundaries[0][i] + " " + boundaries[1][i] + " " + sigma[i]);
        }
        PointValuePair result = calcR.refine(guesses, boundaries[0], boundaries[1], sigma);
        double[] pars = result.getPoint();
        /*
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                System.out.printf(" %3d", map[i][j]);
            }
            System.out.println("");
        }

        System.out.print("Fit pars ");
        for (int i = 0; i < pars.length; i++) {
            System.out.printf(" %.3f", pars[i]);
        }
        System.out.println("");
         */
        double aic = calcR.getAICc(pars);
        double rms = calcR.getRMS(pars);
//        System.out.println("rms " + rms);
        int nGroupPars = calcR.getNGroupPars();
        for (int i = 0; i < guesses.length; i++) {
            sigma[i] /= 2.0;
        }

        String[] parNames = calcR.getParNames();
        double[] errEstimates;
        if (nonParBootStrap) {
            errEstimates = calcR.simBoundsBootstrapStream(pars.clone(), boundaries[0], boundaries[1], sigma);
        } else {
            errEstimates = calcR.simBoundsStream(pars.clone(), boundaries[0], boundaries[1], sigma);

        }
        int nNonGroup = parNames.length - nGroupPars;
        List<int[]> stateList = new ArrayList<>();
        List<String> residueNumbers = new ArrayList<>();
        List<CurveFit> curveFits = new ArrayList<>();
//        System.out.println("ning " + nCurves);
        for (int iCurve = 0; iCurve < nCurves; iCurve++) {
            String stateString = resProps.getStateString(states[iCurve]);
            double[] parArray = new double[parNames.length];
            double[] errArray = new double[parNames.length];
            List<ParValueInterface> parValues = new ArrayList<>();
            for (int i = 0; i < nGroupPars; i++) {
                ParValue parValue = new ParValue(resNums[states[iCurve][0]], stateString, parNames[i], pars[i], errEstimates[i]);
                parValues.add(parValue);
                parArray[i] = pars[i];
                errArray[i] = errEstimates[i];
            }
            for (int j = 0; j < nNonGroup; j++) {
                int k = map[iCurve][nGroupPars + j];
                ParValue parValue = new ParValue(resNums[states[iCurve][0]], stateString, parNames[nGroupPars + j], pars[k], errEstimates[k]);
                parValues.add(parValue);
                parArray[nGroupPars + j] = pars[k];
                errArray[nGroupPars + j] = errEstimates[k];
            }
//            System.out.println("res " + resNums[states[iCurve][0]] + " " + parValues.toString());
            stateList.add(states[iCurve].clone());
            residueNumbers.add(resNums[states[iCurve][0]]);

            HashMap<String, Double> parMap = new HashMap<>();
            for (ParValueInterface parValue : parValues) {
                parMap.put(parValue.getName(), parValue.getValue());
                parMap.put(parValue.getName() + ".sd", parValue.getError());
            }
            parMap.put("AIC", aic);
            parMap.put("RMS", rms);
            parMap.put("Equation", 1.0 + equationNameList.indexOf(eqn));

            PlotEquation plotEquation = new PlotEquation(eqn, parArray, errArray, usedFields);
            CurveFit curveFit = new CurveFit(stateString, resNums[states[iCurve][0]], parMap, plotEquation);
            curveFits.add(curveFit);
        }
        CPMGFitResult fitResult = new CPMGFitResult(parNames, curveFits, eqn, nGroupPars, aic, rms);
        return fitResult;
    }

}