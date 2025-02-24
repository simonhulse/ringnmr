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

import static java.lang.Math.abs;

import java.io.InputStream;
import java.io.IOException;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

import static org.comdnmr.modelfree.RelaxEquations.GAMMA_MAP;
import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.util.DataUtil;
import org.comdnmr.util.Utilities;
import org.comdnmr.util.traindata.DataGenerator;

import org.tensorflow.GraphOperation;
import org.tensorflow.SavedModelBundle;
import org.tensorflow.Tensor;
import org.tensorflow.ndarray.FloatNdArray;
import org.tensorflow.ndarray.NdArrays;
import org.tensorflow.ndarray.Shape;
import org.tensorflow.types.TFloat32;

import static org.comdnmr.modelfree.RelaxEquations.GAMMA_MAP;

import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.util.DataUtil;
import org.comdnmr.util.Utilities;
import static org.comdnmr.util.UnzipUtil.unzip;
import static org.comdnmr.util.traindata.DataGenerator.getInterpolatedProfile;

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
        public double[] updateGuess(int i, double[] current, double[] addition, int[] map) {
            // kex is always at index 0
            double kexCurrent = current[0];
            double kexResidue = addition[0];
            double kex = i * kexCurrent / (i + 1) + kexResidue / (i + 1);
            current[0] = kex;

            // R2s
            for (int idx = 1; idx < map.length - 1; idx++) {
                current[map[idx]] = addition[idx];
            }

            // deltaPPMmin
            current[map[map.length - 1]] = addition[map.length - 1];

            return current;
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

        @Override
        public double[] updateGuess(int i, double[] current, double[] addition, int[] map) {
            // kex is always at index 0
            double kexCurrent = current[0];
            double kexResidue = addition[0];
            double kex = i * kexCurrent / (i + 1) + kexResidue / (i + 1);
            current[0] = kex;

            // pa is always at index 1
            double paCurrent = current[1];
            double paResidue = addition[1];
            double pa = i * paCurrent / (i + 1) + paResidue / (i + 1);
            current[1] = pa;

            // R2s
            for (int idx = 2; idx < map.length - 2; idx++) {
                current[map[idx]] = addition[idx];
            }

            // deltaX
            current[map[map.length - 2]] = addition[map.length - 2];

            // deltaH
            current[map[map.length - 1]] = addition[map.length - 1];

            return current;
        }

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
        public double[] updateGuess(int i, double[] current, double[] addition, int[] map) {
            // kex is always at index 0
            double kexCurrent = current[0];
            double kexResidue = addition[0];
            double kex = i * kexCurrent / (i + 1) + kexResidue / (i + 1);
            current[0] = kex;

            // pa is always at index 1
            double paCurrent = current[1];
            double paResidue = addition[1];
            double pa = i * paCurrent / (i + 1) + paResidue / (i + 1);
            current[1] = pa;

            // R2s
            for (int idx = 2; idx < map.length - 1; idx++) {
                current[map[idx]] = addition[idx];
            }

            // deltaX
            current[map[map.length - 1]] = addition[map.length - 1];

            return current;
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

    // >>> Neural Network stuff >>>

    // Placeholders are for:
    // 1. Enum name (CPMGFAST, CPMGSLOW, CPMGMQ)
    // 2. Number of separate profiles (1, 2, 3)
    final String networkPathTemplate = "CPMGEquation/%s/%d/";
    // Largest permitted profile size.
    // Currently hard-coded, though might be ideal to do this in the future:
    //
    // Determined at runtime by analysing the size of the input layer of a
    // neural network.
    // See `fetchMaxProfileSize()`
    //
    final int maxProfileSize = 16;
    // nuCPMG values to interpolate of using neural networks which rely on
    // interpolation
    final List<Double> networkInterpolationXs = Arrays.asList(
        8.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0,
        100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0);

    // <<< Neural Network stuff <<<

    CPMGEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
        // TODO: see `fetchMaxProfileSize()` below
        // this.maxProfileSize = fetchMaxProfileSize();
    }

    // TODO: may be worth looking at in the future
    // private Integer fetchMaxProfileSize() {
    //     try {
    //         // Setup network loader
    //         NetworkLoader networkLoader = NetworkLoader.getNetworkLoader();
    //         // Open network which considers a single CPMG profile
    //         String networkPath = getNetworkPath("CPMGFAST", 1);
    //         SavedModelBundle network = networkLoader.fetchNetwork(networkPath);

    //         Iterator<GraphOperation> ops = network.graph().operations();
    //         while (ops.hasNext()) {
    //             GraphOperation op = ops.next();
    //             System.out.println(String.format("op.name(): %s", op.name()));
    //         }
    //         // // FIXME: `Shape.size(int i)` has been deprecated.
    //         // // This will likely not work in the future.
    //         // int inputSize = (int) network.graph().operation("serving_default").output(0).shape().size(1);

    //         // For a single profile network: inputSize = 2 * maxProfileSize + 2 + 1
    //         // (2 * maxProfileSize) for x- and y-data
    //         // 2 For fieldH and fieldX
    //         // 1 For tau
    //         return 16;
    //         // return (inputSize - 3) / 2;
    //     } catch (IOException exe) {
    //         return null;
    //     }
    // }

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

    private String getNetworkPath(int nProfiles) {
        return String.format(networkPathTemplate, getName(), nProfiles);
    }

    private String getNetworkPath(String name, int nProfiles) {
        return String.format(networkPathTemplate, name, nProfiles);
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
            try {
                guess = guessNeuralNetwork(xValues, yValues, map, idNums, nID);
            } catch (IOException exe) {
                // Fallback to guessing uing rubric
                guess = guessRubric(xValues, yValues, map, idNums, nID);
            };
        } else {
            guess = guessRubric(xValues, yValues, map, idNums, nID);
        }
        return guess;
    }

    // Will not be called: Is overridden by all enums (see above)
    public double[] guessRubric(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        return new double[0];
    }

    // Will never be called: Is overridden by all enums (see above)
    public double[] updateGuess(int i, double[] current, double[] addition, int[] map) {
        return current;
    }

    public double[] guessNeuralNetwork(
        double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID
    ) throws IOException {
        NNDataCPMG data = new NNDataCPMG(xValues, yValues, map, idNums, maxProfileSize);
        NetworkLoader networkLoader = NetworkLoader.getNetworkLoader();
        double[] guess = new double[getNPars(map)];

        for (int residue = 0; residue < data.getNResidues(); residue++) {
            TFloat32 input = data.getInput(residue);
            int nProfiles = data.getNProfiles(residue);
            String networkPath = getNetworkPath(nProfiles);
            SavedModelBundle network = networkLoader.fetchNetwork(networkPath);
            TFloat32 residueGuessTF32 = (TFloat32) network.function("serving_default").call(input);

            // TFloat32 -> double[]
            int size = (int) residueGuessTF32.size();
            double[] residueGuess = new double[size];
            for (int i = 0; i < residueGuessTF32.size(); i++) {
                residueGuess[i] = residueGuessTF32.getFloat(0, i);
            }

            int[] residueMap = data.getResidueMap(residue);
            guess = updateGuess(residue, guess, residueGuess, residueMap);
        }
        for (var v : guess) {
            System.out.println(v);
        }

        return guess;
    }

    // TODO: update
    //
    // TFloat32 constructNeuralNetworkInputInterpolate(double[][] xValues, double[] yValues) {
    //     Map<Double, Map<String, List<Double>>> datasets = separateDatasets(xValues, yValues);
    //     Map<Double, List<Double>> interpolatedDatasets = interpolateDatasets(datasets);

    //     // Assuming tau is identical across samples
    //     double tau = xValues[3][0];

    //     int nProfiles = interpolatedDatasets.size();
    //     int nValuesPerProfile = networkInterpolationXs.size();
    //     int inputSize = nProfiles * (nValuesPerProfile + 1) + 1;

    //         field = pInfo.getKey();
    //         profile = dataset.getValue();
    //     float[] inputArray = new float[inputSize];
    //     int idx = 0;

    //     double field;
    //     double value;
    //     List<Double> profile;
    //     for (Map.Entry<Double, List<Double>> dataset : interpolatedDatasets.entrySet()) {
    //         // The Map used is a TreeMap, so the iteration will run with the fields in
    //         // order, as desired
    //         field = dataset.getKey();
    //         profile = dataset.getValue();
    //         for (int i = 0; i < nValuesPerProfile; i++) {
    //             value = profile.get(i);
    //             inputArray[idx++] = (float) value;
    //         }
    //         inputArray[idx++] = (float) field;
    //     }
    //     inputArray[idx] = (float) tau;

    //     FloatNdArray inputNdArray = NdArrays.ofFloats(Shape.of(1, inputSize));
    //     inputNdArray.set(TFloat32.vectorOf(inputArray), 0);
    //     return TFloat32.tensorOf(inputNdArray);
    // }

    // TODO: This should be generic across CPMG/CEST/R1rho
    // TODO: Update
    Map<Double, List<Double>> interpolateDatasets(Map<Double, Map<String, List<Double>>> datasets) {
        Map<Double, List<Double>> result = new TreeMap<>();

        double variable;
        Map<String, List<Double>> xyMap;
        List<Double> xs;
        List<Double> ys;
        List<Double> ysInterpolated;

        for (Map.Entry<Double, Map<String, List<Double>>> dataset : datasets.entrySet()) {
            variable = dataset.getKey();
            xyMap = dataset.getValue();
            xs = xyMap.get("x");
            ys = xyMap.get("y");

            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // TODO: Hack to get this to work for datasets with nuCPMGs that do
            // not extend out to spline interpolation limits
            int xIdx = xs.size() - 1;
            int xInterpIdx = networkInterpolationXs.size() - 1;
            if (xs.get(xIdx) < networkInterpolationXs.get(xInterpIdx)) {
                xs.add(networkInterpolationXs.get(xInterpIdx));
                ys.add(ys.get(xIdx));
            }
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            ysInterpolated = getInterpolatedProfile(ys, xs, networkInterpolationXs);
            result.put(variable, ysInterpolated);
        }

        return result;
    }


    // TODO: might as well remove this method?
    // I wrote this when each NN was for a specific nucleus.
    // Now individual NNs can cope with both 13C and 15N, so this is not
    // actually used anywhere in RING
    private String determineNucleus(double fieldH, double fieldX, double epsilon) {
        double targetGyroRatio = fieldX / fieldH;
        double gyroH = GAMMA_MAP.get("H");

        String element;
        double gyroX;
        double gyroRatio;
        for (Map.Entry<String, Double> info : GAMMA_MAP.entrySet()) {
            element = info.getKey();
            gyroX = info.getValue();
            gyroRatio = abs(gyroX / gyroH);
            if (abs(gyroRatio - targetGyroRatio) < epsilon) {
                if (element == "C") {
                    return "13C";
                }
                else if (element == "N") {
                    return "15N";
                }
            }
        }

        throw new IllegalArgumentException(
            "Could not determine heteronucleus for neural network guesser.");
    }
}


class ProfileInfo implements Comparable<ProfileInfo> {
    double fieldH;
    double fieldX;
    List<Double> nuCPMGs;
    List<Double> profile;

    ProfileInfo(double fieldH, double fieldX) {
        this.fieldH = fieldH;
        this.fieldX = fieldX;
        nuCPMGs = new ArrayList<>();
        profile = new ArrayList<>();
    }

    public void addSample(double nuCPMG, double R2) {
        int idx = 0;
        while (idx < nuCPMGs.size() - 1 && nuCPMG > nuCPMGs.get(idx)) {
            idx++;
        }
        nuCPMGs.add(idx, nuCPMG);
        profile.add(idx, R2);
    }

    public double getFieldH() {
        return fieldH;
    }

    public double getFieldX() {
        return fieldX;
    }

    public double getNuCPMG(int idx) {
        if (idx < size()) return nuCPMGs.get(idx);
        else return 0.0;
    }

    public List<Double> getNuCPMG() {
        return nuCPMGs;
    }

    public double getR2(int idx) {
        if (idx < size()) return profile.get(idx);
        else return 0.0;
    }

    public List<Double> getR2() {
        return profile;
    }

    public int size() {
        return nuCPMGs.size();
    }

    @Override
    public int compareTo(ProfileInfo other) {
        return Double.compare(getFieldH(), other.getFieldH());
    }

    @Override
    public String toString() {
        return String.format(
            "ProfileInfo:\n\tfieldH: %s\n\tfieldX: %s\n\tnuCPMGS: %s\n\tprofile: %s\n",
            getFieldH(), getFieldX(), nuCPMGs, profile);
    }
}

class NNDataCPMG {
    private double[][] xData;
    private double[] yData;
    private int[][] profileMap;
    private int[] idNums;
    private int maxProfileSize;
    private float tau;

    private Map<Integer, Pair<Integer, Integer>> idToProfileMap;
    private List<List<ProfileInfo>> residues;
    private List<TFloat32> inputs;
    private int[][] residueMap;

    public NNDataCPMG(double[][] xData, double[] yData, int[][] profileMap, int[] idNums, int maxProfileSize) {
        this.xData = xData;
        this.yData = yData;
        this.profileMap = profileMap;
        this.idNums = idNums;
        this.maxProfileSize = maxProfileSize;
        this.tau = (float) xData[3][0];

        idToProfileMap = getIdToProfileMap();
        residues = separateDatasets();
        inputs = constructInputs();
        residueMap = constructResidueMap();
    }

    public int getNResidues() {
        return inputs.size();
    }

    public int getNProfiles(int residue) {
        return residues.get(residue).size();
    }

    public TFloat32 getInput(int residue) {
        return inputs.get(residue);
    }

    public int[] getResidueMap(int residue) {
        return residueMap[residue];
    }

    private Map<Integer, Pair<Integer, Integer>> getIdToProfileMap() {
        Map<Integer, Integer> residueTracker = new HashMap<>();
        Map<Integer, Integer> profileTracker = new HashMap<>();
        Map<Integer, Pair<Integer, Integer>> idToProfileMap = new HashMap<>();

        int shiftIdx = profileMap[0].length - 1;
        int residueCount = 0;

        int id;
        int idx1;
        int idx2;

        for (int i = 0; i < profileMap.length; i++) {
            id = profileMap[i][shiftIdx];

            if (!residueTracker.containsKey(id)) {
                residueTracker.put(id, residueCount);
                residueCount += 1;
            }
            idx1 = residueTracker.get(id);

            if (!profileTracker.containsKey(id)) {
                profileTracker.put(id, 0);
            }
            else {
                profileTracker.put(id, profileTracker.get(id) + 1);
            }
            idx2 = profileTracker.get(id);

            idToProfileMap.put(i, Pair.of(idx1, idx2));
        }

        return idToProfileMap;
    }

    private List<List<ProfileInfo>> separateDatasets() {
        List<List<ProfileInfo>> residues = initProfiles();

        ProfileInfo pInfo;
        double fieldH;
        double fieldX;
        double nuCPMG;
        double R2;

        int id;
        Pair<Integer, Integer> indices;
        int residueIdx;
        int profileIdx;

        for (int i = 0; i < xData[0].length; i++) {
            fieldH = xData[2][i];
            fieldX = xData[1][i];
            nuCPMG = xData[0][i];
            R2 = yData[i];

            id = idNums[i];
            indices = idToProfileMap.get(id);
            residueIdx = indices.getLeft();
            profileIdx = indices.getRight();

            ProfileInfo profileInfo = residues.get(residueIdx).get(profileIdx);
            if (profileInfo == null) {
                profileInfo = new ProfileInfo(fieldH, fieldX);
            }

            profileInfo.addSample(nuCPMG, R2);
            residues.get(residueIdx).set(profileIdx, profileInfo);
        }

        return residues;
    }

    private List<List<ProfileInfo>> initProfiles()
    {
        int nResidues = 0;
        int nProfiles = 0;
        int residueIdx = 0;
        int profileIdx = 0;

        for (Pair<Integer, Integer> indices : idToProfileMap.values()) {
            residueIdx = indices.getLeft();
            profileIdx = indices.getRight();

            if (residueIdx >= nResidues) {
                nResidues = residueIdx + 1;
            }
            if (profileIdx >= nProfiles) {
                nProfiles = profileIdx + 1;
            }

        }

        List<List<ProfileInfo>> profiles = new ArrayList<>();
        for (int i = 0; i < nResidues; i++) {
            profiles.add(new ArrayList<>());
            for (int j = 0; j < nProfiles; j++) {
                profiles.get(i).add(null);
            }
        }

        return profiles;
    }

    private List<ProfileInfo> filterXValuesNotInAllProfiles(List<ProfileInfo> dataset) {
        Set<Double> allNuCPMGs = new HashSet<>();
        for (ProfileInfo profile : dataset) {
            allNuCPMGs.addAll(profile.getNuCPMG());
        }

        for (double nu : allNuCPMGs) {
            // Check if in all profiles
            boolean inAllProfiles = true;
            for (ProfileInfo pInfo : dataset) {
                if (!pInfo.getNuCPMG().contains(nu)) {
                    inAllProfiles = false;
                    break;
                }
            }
            if (!inAllProfiles) {
                for (ProfileInfo pInfo : dataset) {
                    List<Double> nuCPMG = pInfo.getNuCPMG();
                    List<Double> profile = pInfo.getR2();
                    int idx = nuCPMG.indexOf(nu);
                    while (idx != -1) {
                        nuCPMG.remove(idx);
                        profile.remove(idx);
                        idx = nuCPMG.indexOf(nu);
                    }
                }
            }
        }
        return dataset;
    }

    private List<TFloat32> constructInputs()
    {
        // Assuming tau is identical across samples

        float[] inputArray;
        FloatNdArray inputNdArray;
        TFloat32 inputTensor;

        List<TFloat32> residueTensors = new ArrayList<>();
        for (List<ProfileInfo> residue : residues) {
            residue = filterXValuesNotInAllProfiles(residue);
            inputArray = buildInputArray(residue);

            inputNdArray = NdArrays.ofFloats(Shape.of(1, inputArray.length));
            inputNdArray.set(TFloat32.vectorOf(inputArray), 0);
            inputTensor = TFloat32.tensorOf(inputNdArray);

            residueTensors.add(inputTensor);
        }

        return residueTensors;
    }

    private float[] buildInputArray(List<ProfileInfo> profiles) {
        int nProfiles = profiles.size();
        int inputSize = nProfiles * 2 * (1 * maxProfileSize + 1) + 1;

        float[] inputArray = new float[inputSize];
        int idx = 0;
        for (ProfileInfo pInfo : profiles) {
            for (int i = 0; i < maxProfileSize; i++) {
                inputArray[idx++] = (float) pInfo.getNuCPMG(i);
            }
            for (int i = 0; i < maxProfileSize; i++) {
                inputArray[idx++] = (float) pInfo.getR2(i);
            }
            inputArray[idx++] = (float) pInfo.getFieldH();
            inputArray[idx++] = (float) pInfo.getFieldX();
        }
        inputArray[idx] = tau;

        return inputArray;
    }

    private int[][] constructResidueMap() {
        Map<Integer, SortedSet<Integer>> lookup = new HashMap<>();
        SortedSet<Integer> residueMapSet;
        for (Map.Entry<Integer, Pair<Integer, Integer>> entry : idToProfileMap.entrySet()) {
            int id = entry.getKey();
            int[] pMap = profileMap[id];

            int residueId = entry.getValue().getLeft();

            if (lookup.containsKey(residueId)) {
                residueMapSet = lookup.get(residueId);
            }
            else {
                // Use of linked hash set to ensure correct iteration ordering
                // First into set = first in iteration
                residueMapSet = new TreeSet<>();
                lookup.put(residueId, residueMapSet);
            }

            for (int i = 0; i < pMap.length; i++) {
                residueMapSet.add(pMap[i]);
            }
        }

        int nResidues = lookup.size();
        int nParams = lookup.get(0).size();
        int[][] residueMap = new int[nResidues][nParams];
        for (int r = 0; r < nResidues; r++) {
            int p = 0;
            for (int mapValue : lookup.get(r)) {
                residueMap[r][p++] = mapValue;
            }
        }
        return residueMap;
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
