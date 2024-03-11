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

import org.comdnmr.util.CoMDPreferences;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Bruce Johnson
 */
public interface CESTEquationType extends EquationType {

    @Override
    default double calculate(double[] par, int[] map, double[] X, int idNum) {
        double[][] x = new double[4][1];
        x[0][0] = X[0];
        x[1][0] = X[1];
        x[2][0] = X[2];
        x[3][0] = X[3];
        double[] y = calculate(par, map, x, idNum);
        return y[0];
    }

    @Override
    default double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        double[] guesses;
        if (Boolean.TRUE.equals(CoMDPreferences.getNeuralNetworkGuess())) {
            guesses = ANNguess(xValues, yValues, map, idNums, nID);
        } else {
            guesses = oldGuess(xValues, yValues, map, idNums, nID);
        }
        return guesses;
    }
    
    default double[] ANNguess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        int nPars = CESTFitFunction.getNPars(map);
        double[] guesses = new double[nPars];
        double[] annGuess = new double[map[0].length];
        
        // ANN GUESSER
        for (int id = 0; id < map.length; id++) {
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            double b1Field = xValues[1][0];
            double tEx = xValues[2][0];
            double field = xValues[3][0];
            List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy, "cest");
            try {
                annGuess = CESTEquations.cestANNGuess(xy, peaks, tEx, b1Field);
            } catch (Exception ex) {
                Logger.getLogger(CESTEquation.class.getName()).log(Level.SEVERE, null, ex);
            }
            for (int j = 0; j < map1.length; j++) {
                guesses[map1[j]] = annGuess[j];
            }
        }
        return guesses;
    }
    
    default double[] oldGuess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        int nPars = CESTFitFunction.getNPars(map);
        double[] guesses = new double[nPars];
        
        for (int id = 0; id < map.length; id++) {
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            int yIndex = xy.length - 1;

            List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy, "cest");
            if (!peaks.isEmpty()) {
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(xy[yIndex], tex, "cest");
                double[][] r2 = CESTEquations.cestR2Guess(peaks, xy[yIndex], "cest");
                guesses[map1[0]] = CESTEquations.cestKexGuess(peaks, "cest"); //112.0; //kex
                guesses[map1[1]] = CESTEquations.cestPbGuess(peaks, xy[yIndex], "cest"); //0.1; //pb
                guesses[map1[2]] = peaks.get(0).position; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map1[3]] = peaks.get(peaks.size() - 1).position; //400 * 2.0 * Math.PI; //deltaA
                guesses[map1[4]] = r1[0]; //2.4; //R1A
                guesses[map1[5]] = r1[1]; //2.4; //R1B
                guesses[map1[6]] = r2[0][0]; //20.0; //R2A
                guesses[map1[7]] = r2[1][0]; //100.0; //R2B
            } else {
                guesses = null;
                break;
            }
        }
        return guesses;
    }

    @Override
    default double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        double[][] boundaries = new double[2][guesses.length];
        for (int id = 0; id < map.length; id++) {
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy, "cest");
            double dAbound = 0;
            double dBbound = 0;
            double field = xValues[xValues.length -1][0];
            if (peaks.size() > 1) {
                dAbound = (peaks.get(0).getWidths()[1] / field) / 2;
                dBbound = (peaks.get(1).getWidths()[1] / field) / 2;
            } else if (peaks.size() == 1) {
                dAbound = (peaks.get(0).getWidths()[1] / field) / 2;
                dBbound = dAbound;
            }
            double tex = xValues[2][0];
            double r1A = guesses[map1[4]];
            double[] r1BouA = CESTEquations.r1Boundaries(r1A, tex, 0.1);
            double r1B = guesses[map1[5]];
            double[] r1BouB = CESTEquations.r1Boundaries(r1B, tex, 0.1);

            boundaries[0][map1[0]] = 1.0; //kex LB
            boundaries[1][map1[0]] = guesses[map1[0]] * 6; //kex UB
            boundaries[0][map1[1]] = 0.01; //pb LB
            boundaries[1][map1[1]] = 0.25; //pb UB //guesses[1] * 4;
            boundaries[0][map1[2]] = guesses[map1[2]] - dAbound; //deltaA LB
            boundaries[1][map1[2]] = guesses[map1[2]] + dAbound; //deltaA UB
            boundaries[0][map1[3]] = guesses[map1[3]] - dBbound; //deltaB LB
            boundaries[1][map1[3]] = guesses[map1[3]] + dBbound; //deltaB UB
            boundaries[0][map1[4]] = r1BouA[0]; //R1A LB
            boundaries[1][map1[4]] = r1BouA[1]; //R1A UB
            boundaries[0][map1[5]] = r1BouB[0]; //R1B LB
            boundaries[1][map1[5]] = r1BouB[1]; //R1B UB
            boundaries[0][map1[6]] = 2.0; //R2A LB
            boundaries[1][map1[6]] = guesses[map1[6]] * 6; //R2A UB
            boundaries[0][map1[7]] = 2.0; //R2B LB
            boundaries[1][map1[7]] = guesses[map1[7]] * 6; //R2B UB
            if (boundaries[1][map1[7]] < 200.0) {
                boundaries[1][map1[7]] = 200.0;
            }
        }

        return boundaries;
    }

    @Override
    default double getRex(double[] pars, int[] map, double field) {
        return 0.0;
    }

    @Override
    default double getKex(double[] pars) {
        return pars[0];
    }

    @Override
    default double getKex(double[] pars, int id) {
        return pars[0];
    }

    @Override
    default int[][] makeMap(int n) {
        int nP = 8;
        int[][] map = new int[n][nP];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < nP; j++) {
                map[i][j] = nP * i + j;
            }
        }
        return map;
    }

    @Override
    default int[][] makeMap(int n, int m) {
        int[][] map = new int[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                map[i][0] = m * i + j;
            }
        }
        return map;
    }

    @Override
    default int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
        return makeMap(1);
    }
}
