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
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.eqnfit;

import org.comdnmr.BasicFitter;
import org.comdnmr.data.ExperimentSet;

import java.util.*;

import org.comdnmr.fit.FitQuality;
import org.comdnmr.util.CoMDOptions;
import org.nmrfx.chemistry.relax.ResonanceSource;

/**
 * @author Bruce Johnson
 */
public interface EquationFitter extends BasicFitter {

    List<String> getEquationNameList();

    FitFunction getFitModel();

    String getExpType();

    Optional<FitResult> doFit(String eqn, double[] sliderGuesses, CoMDOptions options);


    int[] getStateCount();

    int[][] getStates();

    double[] getSimX(int nPts, double xLB, double xUB);

    default double[][] getFields(List<Double>[] xValues, List<Integer> idNums) {
        Map<Integer, double[]> fieldMap = new HashMap<>();
        int nExtra = xValues.length - 1;
        for (int i = 0; i < xValues[0].size(); i++) {
            double[] v = new double[nExtra];
            for (int j = 0; j < nExtra; j++) {
                v[j] = xValues[j + 1].get(i);
            }
            fieldMap.put(idNums.get(i), v);
        }
        double[][] values = new double[fieldMap.size()][nExtra];
        for (int i = 0; i < fieldMap.size(); i++) {
            values[i] = fieldMap.get(i);
        }
        return values;
    }

    default FitResult getResults(EquationFitter fitter, String eqn, String[] parNames, ResonanceSource[] dynSources, int[][] map, int[][] states,
                                 double[][] allExtras, int nGroupPars, double[] pars, double[] errEstimates, FitQuality fitQuality, double[][] simPars,
                                 boolean hasExchange, CurveFit.CurveFitStats curveStats) {
        int nNonGroup = parNames.length - nGroupPars;
        List<CurveFit> curveFits = new ArrayList<>();
        Map<String, double[]> simsMap = new HashMap<>();
        for (int[] ints : map) {
            for (int j = 0; j < ints.length; j++) {
                String key = parNames[j] + " " + (ints[j] - map[0][j]);
                simsMap.put(key, simPars[ints[j]]);
            }
        }
        simsMap.put("fit", simPars[simPars.length - 1]);
        int nCurves = states.length;
        for (int iCurve = 0; iCurve < nCurves; iCurve++) {
            String stateString = ExperimentSet.getStateString(states[iCurve]);
            double[] parArray = new double[parNames.length];
            double[] errArray = new double[parNames.length];
            List<ParValueInterface> parValues = new ArrayList<>();
            for (int i = 0; i < nGroupPars; i++) {
                ParValue parValue = new ParValue(dynSources[states[iCurve][0]], stateString, parNames[i], pars[i], errEstimates[i]);
                parValues.add(parValue);
                parArray[i] = pars[i];
                errArray[i] = errEstimates[i];
            }
            for (int j = 0; j < nNonGroup; j++) {
                int k = map[iCurve][nGroupPars + j];
                ParValue parValue = new ParValue(dynSources[states[iCurve][0]], stateString, parNames[nGroupPars + j], pars[k], errEstimates[k]);
                parValues.add(parValue);
                parArray[nGroupPars + j] = pars[k];
                errArray[nGroupPars + j] = errEstimates[k];
            }

            HashMap<String, Double> parMap = new HashMap<>();
            for (ParValueInterface parValue : parValues) {
                parMap.put(parValue.getName(), parValue.getValue());
                parMap.put(parValue.getName() + ".sd", parValue.getError());
            }
            parMap.put("AIC", fitQuality.aic());
            parMap.put("AICc", fitQuality.aicc());
            parMap.put("RMS", fitQuality.rms());
            parMap.put("rChiSq", fitQuality.rChiSq());
            FitFunction model = getFitModel();

            parMap.put("Equation", 1.0 + fitter.getEquationNameList().indexOf(eqn));

            double[] extras = allExtras == null ? new double[0] : allExtras[iCurve];

            PlotEquation plotEquation = new PlotEquation(fitter.getExpType(), eqn, parArray, errArray, extras);
            CurveFit curveFit = new CurveFit(stateString, dynSources[states[iCurve][0]], parMap, plotEquation);
            curveFits.add(curveFit);
        }
        return new FitResult(parNames, curveFits, eqn, nGroupPars, fitQuality, simsMap, hasExchange, curveStats);
    }
}
