//ringnmr/src/main/java/org/comdnmr/eqnfit/CPMGEquation.java
//Simon Hulse
//simonhulse@protonmail.com
//Last Edited: Wed 16 Oct 2024 11:16:34 AM EDT

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

import org.apache.commons.lang3.ArrayUtils;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.util.DataUtil;
import org.comdnmr.util.Utilities;
import org.comdnmr.util.traindata.DataGenerator;

import org.tensorflow.SavedModelBundle;
import org.tensorflow.Tensor;
import org.tensorflow.ndarray.FloatNdArray;
import org.tensorflow.ndarray.NdArrays;
import org.tensorflow.ndarray.Shape;
import org.tensorflow.types.TFloat32;

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Bruce Johnson
 */
public enum CPMGEquation implements EquationType {

    NOEX("NOEX", 0, "R2") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            return par[map[0]];
        }

        @Override
        public double[] guessRubric(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[] guesses = new double[getNPars(map)];
            for (int id = 0; id < map.length; id++) {
                double mean = DataUtil.getMeanValue(yValues, idNums, id);
                guesses[map[id][0]] = mean;
            }
            return guesses;
        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[][] boundaries = new double[2][guesses.length];
            for (int id = 0; id < guesses.length; id++) {
                boundaries[0][id] = 0.0;
                boundaries[1][id] = guesses[id] * 4;
            }
            return boundaries;
        }

        @Override
        public double getKex(double[] pars) {
            return 0.0;
        }

        @Override
        public double getKex(double[] pars, int id) {
            return getKex(pars);
        }

        @Override
        public double getRex(double[] pars, int[] map, double field) {
            return 0.0;
        }

        @Override
        public int[][] makeMap(int n) {
            int[][] map = new int[n][1];
            for (int i = 0; i < n; i++) {
                map[i][0] = i;
            }
            return map;
        }

        @Override
        public int[][] makeMap(int n, int m) {
            return makeMap(n);
        }

        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int n = states.length;
            int[][] map = new int[n][1];
            for (int i = 0; i < n; i++) {
                map[i][0] = CPMGFitter.getMapIndex(states[i], stateCount, r2Mask);
            }
            return map;
        }
    },

    CPMGFAST("CPMGFAST", 1, "Kex", "R2", "dPPMmin") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            double kEx = par[map[0]];
            double R2 = par[map[1]];
            double dPPMmin = par[map[2]];
            double vu = x[0];
            double field = x[1];
            double value;
            if (kEx <= 0.0) {
                value = R2;
            } else {
                double tauCP = 1.0 / (2.0 * vu);
                double dPPMMinRad = 2.0 * Math.PI * dPPMmin * field;
                double Rex = dPPMMinRad * dPPMMinRad / 4.0 / kEx;
                value = R2 + Rex * (1 - 2.0 * FastMath.tanh(0.5 * kEx * tauCP) / (kEx * tauCP));
            }
            return value;
        }

        @Override
        public double[] guessRubric(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[] guesses = new double[getNPars(map)];
            double kExSum = 0.0;
            for (int id = 0; id < map.length; id++) {
                double minY = DataUtil.getMinValue(yValues, idNums, id);
                double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                double r2 = minY * 0.95;
                double rex = maxY - minY;
                if (rex < 0.0) {
                    rex = 0.0;
                }
                double field = xValues[1][id];
                guesses[map[id][1]] = r2;
                double tauMid = 1.0 / (2.0 * vMid);
                double kEx = 1.915 / (0.5 * tauMid);
                double dPPMMinRad = Math.sqrt(4.0 * rex / (field * field) * kEx);
                double dPPMMin = dPPMMinRad / (2.0 * Math.PI);
                guesses[map[id][2]] = dPPMMin;
                if (rex >= 0) {
                    kExSum += kEx; // 1.915 comes from solving equation iteratively at tcp rex 0.5 half max
                }
            }
            guesses[0] = kExSum / map.length;
            if (guesses[0] > CoMDPreferences.getCPMGMaxFreq()) {
                guesses[0] = CoMDPreferences.getCPMGMaxFreq() * 0.9;
            }
        return guesses;
        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[][] boundaries = new double[2][guesses.length];
            for (int[] ints : map) {
                int iPar = ints[0];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = Math.min(guesses[iPar] * 4, CoMDPreferences.getCPMGMaxFreq());
                iPar = ints[1];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
                iPar = ints[2];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
            }
            return boundaries;
        }

        @Override
        public double getRex(double[] pars, int[] map, double field) {
            double dPPMmin = pars[map[2]];

            double kEx = pars[0];
            double dPPMMinRad = 2.0 * Math.PI * dPPMmin * field;
            return dPPMMinRad * dPPMMinRad / 4.0 / kEx;
        }

        @Override
        public double getKex(double[] pars) {
            return pars[0];
        }

        @Override
        public double getKex(double[] pars, int id) {
            return getKex(pars);
        }

        @Override
        public int[][] makeMap(int n) {
            int[][] map = new int[n][3];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 2 * i + 1;
                map[i][2] = 2 * i + 2;
            }
            return map;
        }

        @Override
        public int[][] makeMap(int n, int m) {
            return makeMap(n);
        }

        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int n = states.length;
            int[][] map = new int[n][3];
            int lastCount = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
            }
            lastCount++;
            int maxIndex = 0;
            for (int i = 0; i < n; i++) {
                map[i][1] = CPMGFitter.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
                maxIndex = Math.max(map[i][1], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][2] = CPMGFitter.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
            }
            return map;
        }
    },

    CPMGMQ("CPMGMQ", 2, "kEx", "pA", "R2", "deltaCPPM", "deltaHPPM") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            // See the following DOI:
            // 10.1021/ja039587i
            // References to equations and comments with LaTeX notation relate to this paper
            double kEx = par[map[0]];
            double pA = par[map[1]];
            double R2 = par[map[2]];
            double deltaCPPM = par[map[3]];
            double deltaHPPM = par[map[4]];

            double pB = 1.0 - pA;
            double vcpmg = x[0];
            double fieldX = x[1];
            double fieldH = x[2];
            double tau = x[3];

            double deltaC = 2.0 * Math.PI * deltaCPPM * fieldX;
            double deltaH = fieldH > 1.0e-6 ?  2.0 * Math.PI * deltaHPPM * fieldH : 0.0;

            // TODO: Need to get number of CPMG cycles (n)
            // TODO: or the total time of the CPMG elemnt (T) => n = T / (4 * delta)
            // N.B. (2 * delta) is the time between successive 13C 180 pulses
            double delta = 1.0 / (4.0 * vcpmg);

            // >>> Building lambda1 (3.2 - 3.6) >>>

            // num1: (p_A - p_B)k_{ex} + i \Delta \omega_H
            Complex num1 = new Complex((pA - pB) * kEx, deltaH);

            // 3.6
            Complex zeta = num1.multiply(-2.0 * deltaC);

            // 3.5
            Complex Psi = num1
                .pow(2.0)
                .subtract(Math.pow(deltaC, 2.0))
                .add(4.0 * pA * pB * Math.pow(kEx, 2.0));

            // num2: \sqrt{\Psi^2 + \zeta^2}
            Complex num2 = Psi
                .pow(2.0)
                .add(zeta.pow(2.0))
                .sqrt();

            // 3.4
            Complex etaPlus = num2
                .add(Psi)
                .sqrt()
                .multiply(Math.sqrt(2.0) * delta);
            Complex etaMinus = num2
                .subtract(Psi)
                .sqrt()
                .multiply(Math.sqrt(2.0) * delta);

            // 3.3
            Complex DPlus = Psi
                .add(2.0 * Math.pow(deltaC, 2.0))
                .divide(num2)
                .add(1.0)
                .multiply(0.5);
            Complex DMinus = Psi
                .add(2.0 * Math.pow(deltaC, 2.0))
                .divide(num2)
                .subtract(1.0)
                .multiply(0.5);

            // num3: D_{+} \cosh \eta_{+} - D_{-} \cos \eta_{-}
            Complex num3 = etaPlus
                .cosh()
                .multiply(DPlus)
                .subtract(etaMinus
                    .cos()
                    .multiply(DMinus)
                );

            // 3.2
            // Using \cosh^{-1}(z) = \ln (z + \sqrt{z + 1}\sqrt{z - 1})
            Complex lambda1 = num3
                .add(num3
                    .add(1.0)
                    .sqrt()
                    .multiply(num3
                        .subtract(1.0)
                        .sqrt()
                    )
                )
                .log()
                .divide(-2 * delta)
                .add(kEx)
                .multiply(0.5)
                .add(R2);

            // <<< Building lambda1 (3.2 - 3.6) <<<

            // >>> Building Q (3.7 - 3.10) >>>

            // 3.10
            Complex dPlus = new Complex(deltaH + deltaC, kEx);
            Complex dMinus = new Complex(deltaH + deltaC, -kEx);
            Complex zPlus = new Complex(deltaH - deltaC, kEx);
            Complex zMinus = new Complex(deltaH - deltaC, -kEx);

            // num4: i k_{ex} \sqrt{p_A p_B}
            Complex num4 = new Complex(0.0, kEx * Math.sqrt(pA * pB));

            // 3.9
            Complex mZ = dMinus
                .subtract(dMinus
                    .multiply(delta)
                    .sin()
                    .divide(dMinus
                        .add(zMinus)
                        .multiply(delta)
                        .sin()
                    )
                    .multiply(2.0 * deltaC)
                )
                .multiply(num4
                    .divide(dMinus
                        .multiply(zMinus)
                    )
                )
                .multiply(-1.0);

            // 3.8
            Complex mD = zPlus
                .add(zPlus
                    .multiply(delta)
                    .sin()
                    .divide(dPlus
                        .add(zPlus)
                        .multiply(delta)
                        .sin()
                    )
                    .multiply(2.0 * deltaC)
                )
                .multiply(num4
                    .divide(dPlus
                        .multiply(zPlus)
                    )
                );

            // 3.7
            double Q = mD
                .pow(2.0)
                .multiply(-1.0)
                .add(1.0)
                .add(mD.multiply(mZ))
                .subtract(mZ.pow(2.0))
                .add(mD
                    .add(mZ)
                    .multiply(0.5 * Math.sqrt(pB / pA))
                )
                .getReal();

            // <<< Building Q (3.7 - 3.10) <<<

            // 3.1
            if (tau > 1.0e-6) {
                return lambda1.getReal() - Math.log(Q) / tau;
            } else {
                return lambda1.getReal();
            }
        }

        // TODO
        @Override
        public double[] guessRubric(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[] guesses = new double[getNPars(map)];
            double kExSum = 0.0;
            double pa = 0.95;
            for (int id = 0; id < map.length; id++) {
                double minY = DataUtil.getMinValue(yValues, idNums, id);
                double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                double r2 = minY * 0.95;
                double rex = maxY - r2;
                double tauMid = 1.0 / (2.0 * vMid);
                double kex = 1.915 / (0.5 * tauMid); // 1.915 comes from solving equation iteratively at tcp rex 0.5 half max
                if (kex > CoMDPreferences.getCPMGMaxFreq()) {
                    kex = CoMDPreferences.getCPMGMaxFreq() * 0.9;
                }
                double fieldX = xValues[1][id];
                double dw2 = rex / (pa * (1.0 - pa)) * kex;
                double dPPMC = Math.sqrt(dw2) / (2.0 * Math.PI) / fieldX;
                double dPPMH = 0.1;
                guesses[map[id][2]] = r2;
                guesses[map[id][3]] = dPPMC;
                guesses[map[id][4]] = dPPMH;
                kExSum += kex;
            }
            guesses[0] = kExSum / map.length;
            guesses[1] = pa;
            return guesses;
        }
        int idxInterp = 0;

        // TODO
        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[][] boundaries = new double[2][guesses.length];
            // "Kex", "pA", "R2", "dPPMH, dPPMC"
            for (int[] ints : map) {
                int iPar = ints[0];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = Math.min(guesses[iPar] * 4, CoMDPreferences.getCPMGMaxFreq());
                iPar = ints[1];
                boundaries[0][iPar] = 0.5;
                boundaries[1][iPar] = 0.999;
                iPar = ints[2];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
                iPar = ints[3];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
                iPar = ints[4];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
            }
            return boundaries;
        }

        // TODO
        @Override
        public double getRex(double[] pars, int[] map, double fields) {
            double[] x = {0, fields, 0.0, 0.0};
            x[0] = 10.0;
            double y0 = calculate(pars, map, x, 0);
            x[0] = 1.0e4;
            double y1 = calculate(pars, map, x, 0);
            return y0 - y1;
        }

        // TODO
        @Override
        public double getKex(double[] pars) {
            return pars[0];
        }

        // TODO
        @Override
        public double getKex(double[] pars, int id) {
            return pars[0];
        }

        // TODO
        @Override
        public int[][] makeMap(int n) {
            int[][] map = new int[n][5];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                map[i][2] = 3 * i + 2;
                map[i][3] = 3 * i + 3;
                map[i][4] = 3 * i + 4;
            }
            return map;
        }

        // TODO
        @Override
        public int[][] makeMap(int n, int m) {
            int[][] map = new int[n][5];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                map[i][2] = 3 * i + 2;
                map[i][3] = 3 * i + 3;
                map[i][4] = 3 * i + 4;
            }
            return map;
        }

        // TODO
        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int n = states.length;
            int[][] map = new int[n][5];
            int lastCount;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
            }
            int maxIndex = 0;
            lastCount = 2;
            for (int i = 0; i < n; i++) {
                map[i][2] = CPMGFitter.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
                maxIndex = Math.max(map[i][2], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][3] = CPMGFitter.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
                maxIndex = Math.max(map[i][3], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][4] = CPMGFitter.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
            }
            return map;
        }
    },

    CPMGSLOW("CPMGSLOW", 2, "Kex", "pA", "R2", "dPPM") {

        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            double nu = x[0];
            double field = x[1];
            double kEx = par[map[0]];
            double pA = par[map[1]]; // p1-p2
            double r2 = par[map[2]];
            double dPPM = par[map[3]];
            double pB = 1.0 - pA;
            double pDelta = pA - pB;
            double dW = dPPM * field * 2.0 * Math.PI;
            double tauCP = 1.0 / (2.0 * nu);
            double psi = (pDelta * kEx) * (pDelta * kEx) - dW * dW + 4.0 * pA * pB * kEx * kEx;
            double zeta = -2.0 * dW * kEx * pDelta;
            double eta1 = Math.sqrt(psi * psi + zeta * zeta);
            double etaP = (1.0 / Math.sqrt(2.0)) * tauCP * Math.sqrt(eta1 + psi);
            double etaM = (1.0 / Math.sqrt(2.0)) * tauCP * Math.sqrt(eta1 - psi);
            double d1 = (psi + 2.0 * dW * dW) / Math.sqrt(psi * psi + zeta * zeta);
            double dP = 0.5 * (d1 + 1);
            double dM = 0.5 * (d1 - 1);
            double ch = dP * Math.cosh(etaP) - dM * Math.cos(etaM);
            double rexContrib = 0.5 * (kEx - (1.0 / tauCP) * FastMath.acosh(ch));
            return r2 + rexContrib;
        }

        @Override
        public double[] guessRubric(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[] guesses = new double[getNPars(map)];
            double kExSum = 0.0;
            double pa = 0.95;
            for (int id = 0; id < map.length; id++) {
                double minY = DataUtil.getMinValue(yValues, idNums, id);
                double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                double r2 = minY * 0.95;
                double rex = maxY - r2;
                double tauMid = 1.0 / (2.0 * vMid);
                double kex = 1.915 / (0.5 * tauMid); // 1.915 comes from solving equation iteratively at tcp rex 0.5 half max
                if (kex > CoMDPreferences.getCPMGMaxFreq()) {
                    kex = CoMDPreferences.getCPMGMaxFreq() * 0.9;
                }
                double field = xValues[1][id];
                double dw2 = rex / (pa * (1.0 - pa)) * kex;
                double dPPM = Math.sqrt(dw2) / (2.0 * Math.PI) / field;
                guesses[map[id][2]] = r2;
                guesses[map[id][3]] = dPPM;
                kExSum += kex;
            }
            guesses[0] = kExSum / map.length;
            guesses[1] = pa;
            return guesses;
        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[][] boundaries = new double[2][guesses.length];
            for (int[] ints : map) {
                int iPar = ints[0];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = Math.min(guesses[iPar] * 4, CoMDPreferences.getCPMGMaxFreq());
                iPar = ints[1];
                boundaries[0][iPar] = 0.5;
                boundaries[1][iPar] = 0.999;
                iPar = ints[2];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
                iPar = ints[3];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
            }
            return boundaries;
        }

        @Override
        public double getRex(double[] pars, int[] map, double field) {
            double[] x = {0.0, field};
            x[0] = 10.0;
            double y0 = calculate(pars, map, x, 0);
            x[0] = 1.0e4;
            double y1 = calculate(pars, map, x, 0);
            return y0 - y1;
        }

        @Override
        public double getKex(double[] pars) {
            return pars[0];
        }

        @Override
        public double getKex(double[] pars, int id) {
            return pars[0];
        }

        @Override
        public int[][] makeMap(int n) {
            int[][] map = new int[n][4];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                map[i][2] = 2 * i + 2;
                map[i][3] = 2 * i + 3;
            }
            return map;
        }

        @Override
        public int[][] makeMap(int n, int m) {
            int[][] map = new int[n][4];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                map[i][2] = 2 * i + 2;
                map[i][3] = 2 * i + 3;
            }
            return map;
        }

        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int n = states.length;
            int[][] map = new int[n][4];
            int lastCount;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
            }
            int maxIndex = 0;
            lastCount = 2;
            for (int i = 0; i < n; i++) {
                map[i][2] = CPMGFitter.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
                maxIndex = Math.max(map[i][2], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][3] = CPMGFitter.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
            }
            return map;
        }
    };

    final String equationName;
    final int nGroupPars;
    final String[] parNames;
    // nuCPMG values used to generate training data for neural network guesser.
    // Required to interpolate data acquired with different nuCPMG values.
    // See Simon's `ringguess-java` repo for details:
    // https://github.com/bjohnsonlab/ringguess-java
    // Specifically, see the specification of `n_cycles` and `tau` in `config.json`.
    final List<Double> interpolationXs = Arrays.asList(10.0, 20.0, 50.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1100.0);
    final String guesserPathTemplate = "/data/ringguess/CPMGEquation-%s/models/model_%d/";

    CPMGEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    String getGuesserPath(int nVariable) {
        return String.format(guesserPathTemplate, getName(), nVariable);
    }

    public static String[] getEquationNames() {
        String[] equationNames = {"NOEX", "CPMGFAST", "CPMGSLOW", "CPMGMQ"};
        return equationNames;
    }

    public String getName() {
        return equationName.toUpperCase();
    }

    public String[] getParNames() {
        return parNames;
    }

    public int getNPars(int[][] map) {
        return CPMGFitFunction.getNPars(map);
    }

    public int getNProfiles(int[] idNums) {
        Set<Integer> idNumsSet = new HashSet<Integer>();
        for (int i = 0; i < idNums.length; i++) {
            idNumsSet.add(idNums[i]);
        }
        return idNumsSet.size();
    }

    public int getNGroupPars() {
        return nGroupPars;
    }

    @Override
    public double getMinX() {
        return 5.0;
    }

    @Override
    public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        double[] guess;
        if (
            Boolean.TRUE.equals(CoMDPreferences.getNeuralNetworkGuess())
            && (getName() != "NOEX")
        ) {
            guess = guessNeuralNetwork(xValues, yValues, map, idNums, nID);
        } else {
            guess = guessRubric(xValues, yValues, map, idNums, nID);
        }
        return guess;
    }



    // Will not be called: Is overridden by all enums of this class (see above)
    public double[] guessRubric(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        return new double[0];
    }

    public double[] guessNeuralNetwork(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        TFloat32 input = constructNeuralNetworkInput(xValues, yValues);
        int nProfiles = getNProfiles(idNums);
        int nPars = getNPars(map);
        double[] guess = runNeuralNetwork(input, nProfiles);
        System.out.println(
            String.format(
                "guess:\n%s",
                Arrays.toString(guess)
            )
        );
        return guess;
    }

    Map<Double, List<List<Double>>> separateDatasets(double[][] xValues, double[] yValues) {
        Map<Double, List<List<Double>>> map = new TreeMap<>();
        for (int i = 0; i < xValues[0].length; i++) {
            Double nuCPMG = xValues[0][i];
            Double fieldH = xValues[2][i];
            Double profile = yValues[i];
            List<Double> pair = Arrays.asList(nuCPMG, profile);
            if (map.containsKey(fieldH)) {
                // Add to list at sorted position.
                // Could use binary search to be more efficient O(log n) vs O(n),
                // but the sizes of the datasets means a noticable performance gain
                // probably wouldn't be realised.
                List<List<Double>> profileValues = map.get(fieldH);
                int idx = 0;
                while (idx < profileValues.size() && nuCPMG > profileValues.get(idx).get(0)) {
                    idx++;
                }
                if (idx == profileValues.size()) {
                    profileValues.add(pair);
                } else {
                    profileValues.add(idx, pair);
                }
            } else {
                List<List<Double>> newList = new ArrayList<>();
                newList.add(Arrays.asList(nuCPMG, profile));
                map.put(fieldH, newList);
            }
        }
        return map;
    }

    // TODO: This should be generic across CPMG/CEST/R1rho
    Map<Double, List<Double>> interpolateDatasets(Map<Double, List<List<Double>>> datasets) {
        Map<Double, List<Double>> result = new TreeMap<>();
        for (Map.Entry<Double, List<List<Double>>> dataset : datasets.entrySet()) {
            double variable = dataset.getKey();
            List<List<Double>> xys = dataset.getValue();
            List<Double> xValues = new ArrayList<>();
            List<Double> yValues = new ArrayList<>();
            for (List<Double> xy : xys) {
                double x = xy.get(0);
                double y = xy.get(1);
                xValues.add(x);
                yValues.add(y);
            }

            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // TODO: Hack to get this to work for datasets with nuCPMGs that do
            // not extend out to spline interpolation requires the
            // interpolating x-values to be within the true x-values
            int xIdx = xValues.size() - 1;
            int xInterpIdx = interpolationXs.size() - 1;
            if (xValues.get(xIdx) < interpolationXs.get(xInterpIdx)) {
                xValues.add(interpolationXs.get(xInterpIdx));
                yValues.add(yValues.get(xIdx));
            }
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            List<Double> yValuesInterpolated = DataGenerator.getInterpolatedProfile(yValues, xValues, interpolationXs);
            result.put(variable, yValuesInterpolated);
        }
        return result;
    }

    TFloat32 constructNeuralNetworkInput(double[][] xValues, double[] yValues) {
        Map<Double, List<List<Double>>> datasets = separateDatasets(xValues, yValues);
        Map<Double, List<Double>> interpolatedDatasets = interpolateDatasets(datasets);
        // Assuming tau is identical across samples
        double tau = xValues[3][0];

        int nProfiles = interpolatedDatasets.size();
        int nValuesPerProfile = interpolatedDatasets.entrySet().iterator().next().getValue().size();
        int inputSize = nProfiles * (nValuesPerProfile + 1) + 1;

        float[] inputArray = new float[inputSize];
        int idx = 0;

        for (Map.Entry<Double, List<Double>> entry : interpolatedDatasets.entrySet()) {
            // The Map used is a TreeMap, so the iteration run with the keys in
            // order
            double field = entry.getKey();
            List<Double> profile = entry.getValue();
            for (int i = 0; i < nValuesPerProfile; i++) {
                double value = profile.get(i);
                inputArray[idx++] = (float) value;
            }
            inputArray[idx++] = (float) field;
        }
        inputArray[idx] = (float) tau;

        FloatNdArray inputNdArray = NdArrays.ofFloats(Shape.of(1, inputSize));
        inputNdArray.set(TFloat32.vectorOf(inputArray), 0);
        return TFloat32.tensorOf(inputNdArray);
    }

    double[] runNeuralNetwork(TFloat32 input, int nVariable) {
        SavedModelBundle network = SavedModelBundle.load(getGuesserPath(nVariable), "serve");
        TFloat32 outputTensor = (TFloat32) network.function("serving_default").call(input);

        int size = (int) outputTensor.size();
        double [] output = new double[size];
        for (int i = 0; i < size; i++) {
            output[i] = outputTensor.getFloat(0, i);
        }

        return output;
    }
}

// TODO ???: Implement additional enum: Ishima-Torchia approximation
//
//        ISHIMA("isima", "R2", "Rex", "PaDw", "Tau") {
//            double calculate(double[] par, double tcp, double field) {
//                /*
//                Ishima and Torchia approximation for skewed populations and all time scales
//                R2(1/tcp)=R2+Rex/(1+Tau*sqrt(144/tcp**4+PaDw**4))
//                 */
//                double R2 = par[0];
//                double Rex = par[1];
//                double PaDw = par[2];
//                double Tau = par[3];
//                double value = R2 + Rex / (1 + Tau * FastMath.sqrt(144.0 / FastMath.pow(tcp, 4) + FastMath.pow(PaDw, 4)));
//                return value;
//            }
//            ...
//        },
