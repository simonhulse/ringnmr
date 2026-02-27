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
import org.comdnmr.util.DataUtil;
import java.util.Arrays;

/**
 *
 * @author mbeckwith
 */
public enum NOEEquation implements EquationType {

    NOE("NOE", 0, "noe") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            double noe = par[map[0]];
//            double R = par[map[1]];
//            double delay = x[0];
            double value = noe; //fixme need better equation
            return value;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            int nPars = ExpFitFunction.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double maxY = DataUtil.getMaxValue(yValues, idNums, id);
//                double vMid = DataUtil.getMidValueZero(yValues, xValues[0], idNums, id);
                guesses[map[id][0]] = maxY;
//                guesses[map[id][1]] = -Math.log(0.5) / vMid;

            }
            return guesses;
        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[][] boundaries = new double[2][guesses.length];
            for (int[] map1 : map) {
                int iPar = map1[0];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
//                iPar = map1[1];
//                boundaries[0][iPar] = 0.0;
//                boundaries[1][iPar] = guesses[iPar] * 4;
            }
            return boundaries;
        }

        @Override
        public double getRex(double[] pars, int[] map, double field) {
            return 0.0;
        }

        @Override
        public double getKex(double[] pars) {
            return 0.0;
        }

        @Override
        public double getKex(double[] pars, int id) {
            return 0.0;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int n = states.length;
            int[][] map = new int[n][2];
//            int lastCount = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
            }
//            lastCount++;
//            int maxIndex = 0;
//            for (int i = 0; i < n; i++) {
//                map[i][1] = NOEFit.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
//                maxIndex = Math.max(map[i][1], maxIndex);
//            }
            return map;
        }
    };
    final String equationName;
    final int nGroupPars;
    String[] parNames;
    double fieldRef;

    @Override
    public String getName() {
        return equationName;
    }

    @Override
    public String[] getParNames() {
        return parNames;
    }

    public void setFieldRef(double[] fields) {
        fieldRef = Arrays.stream(fields).min().getAsDouble();
    }

    public void setFieldRef(double field) {
        fieldRef = field;
    }

    @Override
    public int getNGroupPars() {
        return nGroupPars;
    }

    NOEEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    public static String[] getAllEquationNames() {
        String[] equationNames = {"NOE"};
        return equationNames;
    }

    public static String[] getEquationNames() {
        String[] equationNames = new String[CoMDPreferences.getActiveNOEEquations().size()];
        for (int i = 0; i < equationNames.length; i++) {
            equationNames[i] = CoMDPreferences.getActiveNOEEquations().get(i);
        }
        return equationNames;
    }

    public double getMinX() {
        return 0.0;
    }

}
