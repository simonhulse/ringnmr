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

/**
 *
 * @author Bruce Johnson
 */
public interface EquationType {

    double calculate(double[] par, int[] map, double[] x, int idNum);

    default double[] calculate(double[] par, int[] map, double[][] x, int idNum) {
        double[] yValues = new double[x[0].length];
        double[] ax = new double[x.length];
        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0; j < ax.length; j++) {
                ax[j] = x[j][i];
            }
            yValues[i] = calculate(par, map, ax, idNum);
        }
        return yValues;
    }

    default void constrain(String parName, double[] guesses, double[][] boundaries, int[][] map, int id, double lower, double upper) {
        String[] parNames = getParNames();
        int index = -1;
        for (int i = 0; i < parNames.length; i++) {
            if (parNames[i].equals(parName)) {
                index = i;
                break;
            }
        }
        if (index != -1) {
            int j = map[id][index];
            boundaries[0][j] = lower;
            boundaries[1][j] = upper;
            guesses[j] = (lower + upper) / 2;
        }
    }

    double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID);

    double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID);

    double getRex(double[] pars, int[] map, double field);

    double getKex(double[] pars);

    double getKex(double[] pars, int id);

    default  int[][] makeMap(int nStates) {
        int nPars = getParNames().length;
        int nGroupPars = getNGroupPars();
        int nNonGroupPars = nPars - nGroupPars;
        int[][] map = new int[nStates][nPars];
        for (int i = 0; i < nStates; i++) {
            for (int j=0;j<nGroupPars;j++) {
                map[i][j] = j;
            }
            for (int j=0;j<nNonGroupPars;j++) {
                map[i][nGroupPars + j] = nNonGroupPars * i + nGroupPars + j;
            }
        }
        return map;
    }

    int[][] makeMap(int[] stateCount, int[][] states, int[] mask);

    String[] getParNames();

    int getNGroupPars();

    String getName();

    default double getMinX() {
        return -10000.0;
    }

    default double getMaxX() {
        return 10000.0;
    }

}
