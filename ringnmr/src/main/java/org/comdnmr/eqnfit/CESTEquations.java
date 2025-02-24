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

import org.comdnmr.util.MtxExp;
import org.comdnmr.util.ANNLoader;
import org.comdnmr.util.SavitzkyGolay;
import org.comdnmr.util.Utilities;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import static org.comdnmr.util.Utilities.TWO_PI;
import org.ejml.data.Complex_F64;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.decomposition.eig.WatchedDoubleStepQRDecomposition_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ojalgo.ann.ArtificialNeuralNetwork;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.structure.Access1D;

public class CESTEquations {

    /**
     * CEST no-exchange model.
     *
     * @param X Matrix containing the offset values i.e. CEST irradiation
     * frequency (X[0]), the B1 field values (X[1]), and the Tex values (X[2]).
     * @param deltaA0 deltaA0 value. Ground/Major state peak position.
     * @param R1A R1A value. Longitudinal relaxation rate for Ground/Major
     * state.
     * @param R2A R2A value. Transverse relaxation rate for Ground/Major state.
     * @return CEST intensity array.
     */
    public static double[] noEx(double[][] X, double deltaA0, double R1A, double R2A) {
        double[] omegarf = X[0];
        double[] b1Field = X[1];
        double[] Tex = X[2];
        double[] fields = X[3];

        double trad = Tex[0];
        int size = omegarf.length;

        double[] cos2t = new double[size];
        double[] deltaA = new double[size];
        double[] omegaB1 = new double[size];
        for (int i = 0; i < size; i++) {
            omegaB1[i] = b1Field[i] * 2.0 * Math.PI;
            deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * 2.0 * Math.PI;
            double omegaBar = deltaA[i];
            double we = Math.sqrt(omegaB1[i] * omegaB1[i] + omegaBar * omegaBar);
            cos2t[i] = (omegaBar / we) * (omegaBar / we);
        }

        double[] cest = new double[cos2t.length];

        double[] r1rho = R1RhoEquations.r1rhoPerturbationNoEx(omegaB1, deltaA, R1A, R2A);
        for (int i = 0; i < cos2t.length; i++) {
            cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
        }
        return cest;
    }

    /**
     * CEST exact model. Assumes R1A != R1B and R2A != R2B. Uses matrix
     * exponential. Numerical integration of thermalized Bloch-McConnell rate
     * matrix.
     *
     * @param X Matrix containing the offset values i.e. CEST irradiation
     * frequency (X[0]), the B1 field values (X[1]), and the Tex values (X[2]).
     * @param pb pb value. Population of the Excited/Minor state.
     * @param kex kex value. Rate constant for Major-Minor state exchange.
     * @param deltaA0 deltaA0 value. Ground/Major state peak position.
     * @param deltaB0 deltaB0 value. Excited/Minor state peak position.
     * @param R1A R1A value. Longitudinal relaxation rate for Ground/Major
     * state.
     * @param R1B R1B value. Longitudinal relaxation rate for Excited/Minor
     * state.
     * @param R2A R2A value. Transverse relaxation rate for Ground/Major state.
     * @param R2B R2B value. Transverse relaxation rate for Excited/Minor state.
     * @return CEST intensity array.
     */
    public static double[] exact0(double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R1A, double R1B, double R2A, double R2B) {
        // Performs an exact numerical calculation and returns CEST intensity ratio.
        //
        // X: array containing two arrays:
        //  omegarf: CEST irradiation frequency (ppm)
        //  omega1: B1 field strength (1/s)
        //
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states

        double[] omegarf = X[0];
        double[] b1Field = X[1];
        double[] Tex = X[2];
        double[] fields = X[3];

        // time delay is hard-coded below
        double tdelay = Tex[0];

        double[] cest = new double[omegarf.length];

        double k1 = pb * kex;
        double km1 = (1 - pb) * kex;

        double[] m0 = {0.5, 0, 0, 1 - pb, 0, 0, pb};
        double[] m1 = {0.5, 0, 0, -(1 - pb), 0, 0, -pb};

        double[][] K = {{0, 0, 0, 0, 0, 0, 0},
        {0, -k1, 0, 0, km1, 0, 0},
        {0, 0, -k1, 0, 0, km1, 0},
        {0, 0, 0, -k1, 0, 0, km1},
        {0, k1, 0, 0, -km1, 0, 0},
        {0, 0, k1, 0, 0, -km1, 0},
        {0, 0, 0, k1, 0, 0, -km1}};

        double[][] La = {{0, 0, 0, 0, 0, 0, 0},
        {0, -R2A, 0, 0, 0, 0, 0}, //{0, -R2A, -deltaA, 0, 0, 0, 0}
        {0, 0, -R2A, 0, 0, 0, 0}, //{0, deltaA, -R2A, -omegaB1[i], 0, 0, 0}
        {2 * R1A * (1 - pb), 0, 0, -R1A, 0, 0, 0}, //{2 * R1A * (1 - pb), 0, omegaB1[i], -R1A, 0, 0, 0}
        {0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0}};

        double[][] Lb = {{0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, -R2B, 0, 0}, //{0, 0, 0, 0, -R2B, -deltaB, 0}
        {0, 0, 0, 0, 0, -R2B, 0}, //{0, 0, 0, 0, deltaB, -R2B, -omegaB1[i]}
        {2 * R1B * pb, 0, 0, 0, 0, 0, -R1B}}; //{2 * R1B * pb, 0, 0, 0, 0, omegaB1[i], -R1B}

        DMatrixRMaj La1 = new DMatrixRMaj(La);
        DMatrixRMaj Lb1 = new DMatrixRMaj(Lb);
        DMatrixRMaj K1 = new DMatrixRMaj(K);

        DMatrixRMaj Z = new DMatrixRMaj(La.length, La[0].length);

        for (int i = 0; i < omegarf.length; i++) {
            double omegaB1 = b1Field[i] * TWO_PI;
            double deltaA = (deltaA0 - omegarf[i]) * fields[i] * TWO_PI;
            double deltaB = (deltaB0 - omegarf[i]) * fields[i] * TWO_PI;

            La1.set(1, 2, -deltaA);
            La1.set(2, 1, deltaA);
            La1.set(2, 3, -omegaB1);
            La1.set(3, 2, omegaB1);

            Lb1.set(4, 5, -deltaB);
            Lb1.set(5, 4, deltaB);
            Lb1.set(5, 6, -omegaB1);
            Lb1.set(6, 5, omegaB1);

            CommonOps_DDRM.add(La1, Lb1, Z);
            CommonOps_DDRM.addEquals(Z, K1);

            CommonOps_DDRM.scale(tdelay, Z);

            DMatrixRMaj at = MtxExp.matrixExp(Z);

            double a30 = at.get(3, 0);
            double a33 = at.get(3, 3);
            double a36 = at.get(3, 6);
            double magA = a30 * m0[0] + a33 * m0[3] + a36 * m0[6];
            magA = magA - (a30 * m1[0] + a33 * m1[3] + a36 * m1[6]);
            magA = magA / 2;

            cest[i] = magA;
        }

        return cest;
    }

    /**
     * CEST exact model. Assumes R1A = R1B and R2A != R2B. Uses numerical
     * determination of least negative eigenvalue to calculate CEST.
     *
     * @param X Matrix containing the offset values i.e. CEST irradiation
     * frequency (X[0]), the B1 field values (X[1]), and the Tex values (X[2]).
     * @param pb pb value. Population of the Excited/Minor state.
     * @param kex kex value. Rate constant for Major-Minor state exchange.
     * @param deltaA0 deltaA0 value. Ground/Major state peak position.
     * @param deltaB0 deltaB0 value. Excited/Minor state peak position.
     * @param R1A R1A value. Longitudinal relaxation rate for Ground/Major
     * state.
     * @param R1B R1B value. Longitudinal relaxation rate for Excited/Minor
     * state.
     * @param R2A R2A value. Transverse relaxation rate for Ground/Major state.
     * @param R2B R2B value. Transverse relaxation rate for Excited/Minor state.
     * @return CEST intensity array.
     */
    public static double[] eigenExact1(double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R1A, double R1B, double R2A, double R2B) {
        // Performs an exact numerical calculation and returns CEST intensity ratio.
        // Assumes R1A = R1B.
        //
        // X: array containing two arrays:
        //  omegarf: CEST irradiation frequency (ppm)
        //  omega1: B1 field strength (1/s)
        //
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states

        double[] omegarf = X[0];
        double[] b1Field = X[1];
        double[] Tex = X[2];
        double[] fields = X[3];

        double tdelay = Tex[0];

        double[] cest = new double[omegarf.length];

        double k1 = pb * kex;
        double km1 = (1 - pb) * kex;

        double[][] K = {{-k1, 0, 0, km1, 0, 0},
        {0, -k1, 0, 0, km1, 0},
        {0, 0, -k1, 0, 0, km1},
        {k1, 0, 0, -km1, 0, 0},
        {0, k1, 0, 0, -km1, 0},
        {0, 0, k1, 0, 0, -km1}};

        double[][] La = {{-R2A, 0, 0, 0, 0, 0}, //{-R2A, -deltaA, 0, 0, 0, 0}
        {0, -R2A, 0, 0, 0, 0}, //{deltaA, -R2A, -omegaB1[i], 0, 0, 0}
        {0, 0, -R1A, 0, 0, 0}, //{0, omegaB1[i], -R1A, 0, 0, 0}
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0}};

        double[][] Lb = {{0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, -R2B, 0, 0}, //{0, 0, 0, -R2B, -deltaB, 0}
        {0, 0, 0, 0, -R2B, 0}, //{0, 0, 0, deltaB, -R2B, -omegaB1[i]}
        {0, 0, 0, 0, 0, -R1B}}; //{0, 0, 0, 0, omegaB1[i], -R1B}

        DMatrixRMaj La1 = new DMatrixRMaj(La);
        DMatrixRMaj Lb1 = new DMatrixRMaj(Lb);
        DMatrixRMaj K1 = new DMatrixRMaj(K);

        DMatrixRMaj Z = new DMatrixRMaj(La.length, La[0].length);

        for (int i = 0; i < omegarf.length; i++) {
            double omegaB1 = b1Field[i] * TWO_PI;
            double deltaA = (deltaA0 - omegarf[i]) * fields[i] * TWO_PI;
            double deltaB = (deltaB0 - omegarf[i]) * fields[i] * TWO_PI;
            double omegaBar = (1 - pb) * deltaA + pb * deltaB;
            double we = Math.sqrt(omegaB1 * omegaB1 + omegaBar * omegaBar);

            double cost = omegaBar / we;

            La1.set(0, 1, -deltaA);
            La1.set(1, 0, deltaA);
            La1.set(1, 2, -omegaB1);
            La1.set(2, 1, omegaB1);

            Lb1.set(3, 4, -deltaB);
            Lb1.set(4, 3, deltaB);
            Lb1.set(4, 5, -omegaB1);
            Lb1.set(5, 4, omegaB1);

            CommonOps_DDRM.add(La1, Lb1, Z);
            CommonOps_DDRM.addEquals(Z, K1);

            List<Double> r1rho1 = new ArrayList<>();
            WatchedDoubleStepQRDecomposition_DDRM eig = (WatchedDoubleStepQRDecomposition_DDRM) DecompositionFactory_DDRM.eig(Z.getNumCols(), true, false);
            eig.decompose(Z);
            int numEigVec = eig.getNumberOfEigenvalues();
            for (int v = 0; v < numEigVec; v++) {
                Complex_F64 ZEigVal = eig.getEigenvalue(v);
                if (ZEigVal.isReal()) {
                    r1rho1.add(ZEigVal.getReal());
                }
            }

            double r1rho = Math.abs(Collections.max(r1rho1));

            cest[i] = cost * cost * Math.exp(-tdelay * r1rho);
        }

        return cest;
    }

    /**
     * CEST approximation models. All expect laguerre assume R1A = R1B and R2A
     * != R2B. Laguerre assumes R1A = R1B and R2A = R2B. Laguerre: R1rho
     * approximation to calculate CEST Second-order approximation, assuming
     * R1A=R1B and R2A=R2B Perturbation: Uses Trott perturbation R1rho
     * approximation to calculate CEST Assumes R1A=R1B SD: Uses perturbation
     * R1rho approximation to calculate CEST second-order approximation,
     * assuming R1A=R1B. Averages over B1 inhomogeneity BaldwinKay Uses
     * first-order R1rho approximation to calculate CEST. Assumes R1A=R1B
     *
     * @param X Matrix containing the offset values i.e. CEST irradiation
     * frequency (X[0]), the B1 field values (X[1]), and the Tex values (X[2]).
     * @param pb pb value. Population of the Excited/Minor state.
     * @param kex kex value. Rate constant for Major-Minor state exchange.
     * @param deltaA0 deltaA0 value. Ground/Major state peak position.
     * @param deltaB0 deltaB0 value. Excited/Minor state peak position.
     * @param R1A R1A value. Longitudinal relaxation rate for Ground/Major
     * state.
     * @param R1B R1B value. Longitudinal relaxation rate for Excited/Minor
     * state.
     * @param R2A R2A value. Transverse relaxation rate for Ground/Major state.
     * @param R2B R2B value. Transverse relaxation rate for Excited/Minor state.
     * @return CEST intensity array.
     */
    public static double[] r1rhoApprox(String approx, double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R1A, double R1B, double R2A, double R2B) {

        // X: array containing two arrays:
        //  omegarf: CEST irradiation frequency (ppm)
        //  omega1: B1 field strength (1/s)
        //
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states
        // In the present implementation, the irradiation time is hard-coded below
        double[] omegarf = X[0];
        double[] b1Field = X[1];
        double[] Tex = X[2];
        double[] fields = X[3];

        double trad = Tex[0];

        double pa = 1.0 - pb;
        double[] deltaA = new double[omegarf.length];
        double[] deltaB = new double[omegarf.length];
        double[] omegaBar = new double[omegarf.length];
        double[] omegaB1 = new double[omegarf.length];
        for (int i = 0; i < omegarf.length; i++) {
            omegaB1[i] = b1Field[i] * TWO_PI;
            deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * TWO_PI;
            deltaB[i] = (deltaB0 - omegarf[i]) * fields[i] * TWO_PI;
            omegaBar[i] = pa * deltaA[i] + pb * deltaB[i];
        }

        double[] we = new double[omegaBar.length];
        double[] cos2t = new double[omegaBar.length];
        for (int i = 0; i < omegaBar.length; i++) {
            we[i] = Math.sqrt(omegaB1[i] * omegaB1[i] + omegaBar[i] * omegaBar[i]);
            cos2t[i] = (omegaBar[i] / we[i]) * (omegaBar[i] / we[i]);
        }

        double[] cest = new double[cos2t.length];

        if (null != approx) {
            switch (approx) {
                case "laguerre" -> {
                    double[] r1rho = R1RhoEquations.r1rhoLaguerre(trad, false, omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
                    for (int i = 0; i < cos2t.length; i++) {
                        cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
                    }
                }
                case "trott" -> {
                    double[] r1rho = R1RhoEquations.r1rhoPerturbation(trad, false, omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
                    for (int i = 0; i < cos2t.length; i++) {
                        cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
                    }
                }
                case "trottnoex" -> {
                    double[] r1rho = R1RhoEquations.r1rhoPerturbationNoEx(omegaB1, deltaA, R1A, R2A);
                    for (int i = 0; i < cos2t.length; i++) {
                        cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
                    }
                }
                case "baldwinkay" -> {
                    double[] r1rho = R1RhoEquations.r1rhoBaldwinKay(trad, false, omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
                    for (int i = 0; i < cos2t.length; i++) {
                        cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
                    }
                }
                case "sd" -> {
                    int wlen = 11;
                    double omegaSD = 0.2;   // fractional variation in B1
                    double[] omwt = {0.022, 0.0444, 0.0777, 0.1159, 0.1473, 0.1596, 0.1473, 0.1159, 0.0777, 0.0444, 0.0216};
                    double omwtsum = 0;
                    for (double v : omwt) {
                        omwtsum += v;
                    }
                    for (int i = 0; i < omwt.length; i++) {
                        omwt[i] = omwt[i] / omwtsum;
                    }
                    double[] omegagauss = new double[wlen];
                    for (int i = 0; i < wlen; i++) {
                        omegagauss[i] = -2 * omegaSD + i * (2 * omegaSD - (-2 * omegaSD)) / (wlen - 1);
                    }
                    double[] magA = new double[omegaB1.length];
                    double[] omegatmp = new double[omegaB1.length];
                    for (int i = 0; i < wlen; i++) {
                        for (int j = 0; j < omegatmp.length; j++) {
                            omegatmp[j] = omegaB1[j] * (1 + omegagauss[i]);
                            we[j] = Math.sqrt(omegatmp[j] * omegatmp[j] + omegaBar[j] * omegaBar[j]);
                            cos2t[j] = (omegaBar[j] / we[j]) * (omegaBar[j] / we[j]);
                        }

                        double[] r1rho = R1RhoEquations.r1rhoPerturbation(trad, false, omegatmp, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);

                        for (int j = 0; j < cest.length; j++) {
                            cest[j] = cos2t[j] * Math.exp(-trad * r1rho[j]);
                            magA[j] = magA[j] + omwt[i] * cest[j];
                        }
                    }
                    cest = magA;
                }
                default -> {
                }
            }
        }

        return cest;
    }

    /**
     * Applies a Savitzky-Golay filter to the CEST/R1rho data for smoothing.
     *
     * @param vec Array of CEST/R1rho intensity values to smooth.
     * @param j1 Start index for smoothing in the CEST/R1rho intensity array.
     * @param j2 End index for smoothing in the CEST/R1rho intensity array.
     * @param order Smoothing order. Default is 3.
     * @param smoothSize Size of the smoothing filter. Should be one of 5, 7, 9,
     * ... 25. Depends on intensity array size.
     * @param X Empty array for smoothed values. Same length as intensity array.
     * Returned if an exception is raised.
     * @return Array of smoothed intensity values.
     */
    public static double[] smoothCEST(double[] vec, int j1, int j2, int order, int smoothSize, double[] X) {
        // Applies a Savitzky-Golay filter to the data for smoothing.
        SavitzkyGolay sg = new SavitzkyGolay();
        try {
            return sg.runningSavitzkyGolaySmoothing(vec, j1, j2, order, smoothSize, X);
        } catch (Exception e) {
            System.out.println("Smooth-size should be one of 5,7,9,...,25");
            return X;
        }
    }

    /**
     * Estimates the baseline of the CEST/R1rho data
     *
     * @param vec Array of CEST/R1rho intensity values.
     * @param fitMode String "cest" or "r1rho" specifying CEST or R1rho data.
     * @return Array of the estimate of the baseline in the data.
     */
    public static double[] getBaseline(double[] vec, String fitMode) {
        int winSize = 8;
        double maxValue = (fitMode.equals("r1rho")) ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
        double sDev = 0.0;
        DescriptiveStatistics stat = new DescriptiveStatistics(winSize);
        for (int i = 0; i < vec.length; i++) {
            stat.addValue(vec[i]);
            if (i >= (winSize - 1)) {
                double mean = stat.getMean();
                if ((fitMode.equals("cest") && mean > maxValue) || (fitMode.equals("r1rho") && mean < maxValue)) {
                    maxValue = mean;
                    sDev = stat.getStandardDeviation();
                }
            }
        }
        return new double[]{maxValue, sDev};
    }

    /**
     * Estimates CEST/R1rho peak positions and widths for initial guesses before
     * fitting.
     *
     * @param xyvals Array of the offset values i.e. CEST irradiation frequency.
     * @param fitMode String "cest" or "r1rho" specifying CEST or R1rho data.
     * @return List of CESTPeak objects. The major state peak is listed first,
     * then the minor state peak, if any.
     */
    public static List<CESTPeak> cestPeakGuess(double[][] xyvals, String fitMode) {
        // Estimates CEST peak positions for initial guesses for before fitting.
        double[] xvals = xyvals[0];
        double[] yvals = xyvals[xyvals.length - 1];
        double field = xyvals[xyvals.length-2][0];
        List<CESTPeak> peaks = new ArrayList<>();

        double[] syvals = new double[yvals.length];
        double[] baseValues = getBaseline(yvals, fitMode);
        double baseline = baseValues[0];
        int smoothSize;
        if (yvals.length < 20) {
            smoothSize = 0;
        } else if (yvals.length < 30) {
            smoothSize = 5;
        } else if (yvals.length < 40) {
            smoothSize = 7;
        } else if (yvals.length < 50) {
            smoothSize = 9;
        } else {
            smoothSize = 11;
        }

        if (fitMode.equals("cest") && smoothSize != 0) {
            yvals = smoothCEST(yvals, 0, yvals.length, 3, smoothSize, syvals);
        }

        // A point must have a lower value than this number of points on each
        // side of the point in order to be a peak.
        int nP = 2;
        double baseRatio = (fitMode.equals("r1rho")) ? 1.5 : 3.0;

        // A threshold to use when deciding if a point is deep enough
        // calculated as baseline - a multiple of the standard deviation estimate
        double threshold = baseline - baseValues[1] * baseRatio;
        double yMin = (fitMode.equals("r1rho")) ? Double.MIN_VALUE : Double.MAX_VALUE;

        for (int i = nP; i < yvals.length - nP; i++) {
            if ((fitMode.equals("cest") && yvals[i] < yMin) || (fitMode.equals("r1rho") && yvals[i] > yMin)) {
                yMin = yvals[i];
            }
            if ((fitMode.equals("cest") && yvals[i] < threshold) || (fitMode.equals("r1rho") && yvals[i] > threshold)) {
                boolean ok = true;
                for (int j = i - nP; j <= (i + nP); j++) {
                    if ((fitMode.equals("cest") && yvals[i] > yvals[j]) || (fitMode.equals("r1rho") && yvals[i] < yvals[j])) {
                        ok = false;
                        break;
                    }
                }
                int iCenter = i;
                if (ok) {
                    double halfinten = (baseline - yvals[iCenter]) / 2 + yvals[iCenter];
                    double quarterinten = (baseline - yvals[iCenter]) / 4 + yvals[iCenter];
                    double threequarterinten = 3 * ((baseline - yvals[iCenter]) / 4) + yvals[iCenter];
                    if (fitMode.equals("r1rho")) {
                        halfinten = (yvals[iCenter] - baseline) / 8 + baseline;
                        quarterinten = (yvals[iCenter] - baseline) / 16 + baseline;
                        threequarterinten = 3 * ((yvals[iCenter] - baseline) / 16) + baseline;
                    }

                    double[] widthinten = {halfinten, quarterinten, threequarterinten};
                    double[][] widthpos = new double[widthinten.length][2];
                    boolean[][] foundPos = new boolean[widthinten.length][2];

                    // search from peak center in both directions to find
                    // the peak width50.  Find a value above and below the
                    // half-height and interpolate to get the width50 on each
                    // side.
                    for (int w = 0; w < widthinten.length; w++) {
                        for (int k = 0; k < 2; k++) {
                            int iDir = k * 2 - 1; // make iDir -1, 1
                            int j = iCenter + iDir;
                            double dUp = Double.MAX_VALUE;
                            double dLow = Double.MAX_VALUE;
                            if (fitMode.equals("r1rho")) {
                                dUp = Double.MIN_VALUE;
                                dLow = Double.MIN_VALUE;
                            }
                            int iUp = 0;
                            int iLow = 0;
                            while ((j >= 0) && (j < yvals.length)) {
                                double delta = yvals[j] - widthinten[w];
                                if ((fitMode.equals("cest") && delta < 0.0) || (fitMode.equals("r1rho") && delta > 0.0)) {
                                    if ((fitMode.equals("cest") && Math.abs(delta) < dUp) || (fitMode.equals("r1rho") && Math.abs(delta) > dUp)) {
                                        dUp = Math.abs(delta);
                                        iUp = j;
                                    }
                                } else {
                                    if ((fitMode.equals("cest") && Math.abs(delta) < dLow) || (fitMode.equals("r1rho") && Math.abs(delta) > dLow)) {
                                        dLow = Math.abs(delta);
                                        iLow = j;
                                    }
                                    break;
                                }
                                j += iDir;
                            }
                            foundPos[w][k] = (!fitMode.equals("cest") || ((dLow != Double.MAX_VALUE) && (dUp != Double.MAX_VALUE))) && (!fitMode.equals("r1rho") || ((dLow != Double.MIN_VALUE) && (dUp != Double.MIN_VALUE)));
                            double delta = dLow + dUp;
                            widthpos[w][k] = xvals[iLow] * dUp / delta + xvals[iUp] * dLow / delta;
                        }
                    }
                    double[][] widthtab = new double[widthpos.length][3];
                    if (ok) {
                        double xCenter = xvals[iCenter];
                        double yCenter = yvals[iCenter];
                        for (int w = 0; w < widthpos.length; w++) {
                            if (foundPos[w][0]) {
                                widthtab[w][1] = Math.abs(widthpos[w][0] - xCenter) * field;
                            } else if (foundPos[w][1]) {
                                widthtab[w][1] = Math.abs(widthpos[w][1] - xCenter) * field;
                            } else if (w > 0) {
                                widthtab[w][1] = widthtab[w - 1][1] * 1.3;
                            }
                            if (foundPos[w][1]) {
                                widthtab[w][2] = Math.abs(widthpos[w][1] - xCenter) * field;
                            } else if (foundPos[w][0]) {
                                widthtab[w][2] = Math.abs(widthpos[w][0] - xCenter) * field;
                            } else if (w > 0) {
                                widthtab[w][2] = widthtab[w - 1][2] * 1.3;
                            }
                            widthtab[w][0] = (widthtab[w][1] + widthtab[w][2]);
                        }
                        CESTPeak peak = new CESTPeak(iCenter, xCenter, yCenter, widthtab[0][0], widthtab[0][1], widthtab[0][2], widthtab[1][0], widthtab[1][1], widthtab[1][2], widthtab[2][0], widthtab[2][1], widthtab[2][2], baseline);
                        peaks.add(peak);
                    }
                }
            }
        }

        peaks.sort(Comparator.comparingDouble(CESTPeak::getDepth));

        List<CESTPeak> peaks2 = peaks;
        if (peaks.size() >= 2) {
            peaks2 = peaks.subList(0, 2);
            if (fitMode.equals("r1rho")) {
                peaks2 = peaks.subList(peaks.size() - 2, peaks.size());
            }
        } else if (peaks.size() == 1) {
            // If there is only one peak found add another peak on the side
            // with the largest width50
            peaks2 = peaks.subList(0, 1);
            CESTPeak peak = peaks2.get(0);
            double newCenter;
            if (peak.width50LB > peak.width50UB) {
                newCenter = peak.position + peak.width50LB / field / 2.0;
            } else {
                newCenter = peak.position - peak.width50UB / field / 2.0;
            }
            double newDepth = (baseline + peak.depth) / 2.0;
            CESTPeak newPeak = new CESTPeak(peak.pkInd, newCenter, newDepth, peak.width50, peak.width50LB, peak.width50UB, peak.width25, peak.width25LB, peak.width25UB, peak.width75, peak.width75LB, peak.width75UB, peak.baseline);
            peaks2.add(newPeak);
        }

        peaks2.sort(Comparator.comparingDouble(CESTPeak::getDepth));

        List<CESTPeak> peaks1 = new ArrayList<>();
        if (!peaks2.isEmpty()) {
            double peak1diff = Math.abs(peaks2.get(0).depth - baseline);
            double peak2diff = peak1diff;
            if (peaks2.size() == 2) {
                peak2diff = Math.abs(peaks2.get(1).depth - baseline);
            }
            if (peak1diff < 0.05 || peak2diff < 0.05) {
                int index = 0;
                if (peaks2.get(0).depth > peaks2.get(1).depth && peaks2.get(1).depth != 0) {
                    index = 1;
                }
                peaks1.add(peaks2.get(index));
            } else {
                peaks1.addAll(peaks2);
            }
        } else {
            System.out.println("No peaks found by peak guesser.");
        }
        return peaks1;
    }

    /**
     * Estimates the CEST/R1rho pb values for initial guesses before fitting.
     *
     * @param peaks List of CESTPeak objects from CEST/R1rho peak guesser.
     * @param yvals Array of the CEST intensities.
     * @param fitMode String "cest" or "r1rho" specifying CEST or R1rho data.
     * @return CEST/R1rho pb value estimate.
     */
    public static double cestPbGuess(List<CESTPeak> peaks, double[] yvals, String fitMode) {
        // Estimates CEST pb values from peak intensities for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

        double[] baseValues = getBaseline(yvals, fitMode);
        double baseline = baseValues[0];

        if (peaks.size() > 1) {
            double[] pb = new double[peaks.size() / 2];

            double factor = 4;
            if (fitMode.equals("r1rho")) {
                factor = 40;
            }
            if (peaks.size() == 2) {
                if (peaks.get(0).depth > peaks.get(1).depth) {
                    pb[0] = (Math.abs(baseline - peaks.get(0).depth) / Math.abs(baseline - peaks.get(1).depth)) / factor;
                } else {
                    pb[0] = (Math.abs(baseline - peaks.get(1).depth) / Math.abs(baseline - peaks.get(0).depth)) / factor;
                }
            } else {
                for (int i = 0; i < pb.length; i++) {
                    pb[i] = Math.abs(baseline - peaks.get(2 * i).depth) / Math.abs(baseline - peaks.get(2 * i + 1).depth) / factor;
                }
            }
//            System.out.println("pb guess = " + pb[0]);
            if (pb[0] > 0.25) {
                pb[0] = 0.25;
            }
            return pb[0];
        } else {
            return 0.1;
        }

    }

    /**
     * Estimates the CEST/R1rho R2A and R2B values for initial guesses before
     * fitting.
     *
     * @param peaks List of CESTPeak objects from CEST/R1rho peak guesser.
     * @param yvals Array of the CEST intensities.
     * @param fitMode String "cest" or "r1rho" specifying CEST or R1rho data.
     * @return Matrix of CEST/R1rho R2A and R2B value estimates. R2A is listed
     * first, then R2B.
     */
    public static double[][] cestR2Guess(List<CESTPeak> peaks, double[] yvals, String fitMode) {
        // Estimates CEST R2A and R2B values from peak widths for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

        double pb = cestPbGuess(peaks, yvals, fitMode);

        if (peaks.size() > 1) {
            double[][] r2 = new double[2][peaks.size() / 2];
            double awidth;
            double bwidth;

            double afactor = 1;
            double bfactor = 1;
            if (fitMode.equals("r1rho")) {
                afactor = 12;
                bfactor = 6;
            }
            if (peaks.size() == 2) {
                if (peaks.get(0).depth > peaks.get(1).depth) {
                    awidth = peaks.get(1).width50 / (2 * Math.PI);
                    bwidth = peaks.get(0).width50 / (2 * Math.PI);
                } else {
                    awidth = peaks.get(0).width50 / (2 * Math.PI);
                    bwidth = peaks.get(1).width50 / (2 * Math.PI);
                }
                double kex = (awidth + bwidth) / 2;
                double kb = pb * kex;
                double ka = (1 - pb) * kex;
                r2[0][0] = Math.abs(awidth - ka) / afactor; //R2A
                r2[1][0] = Math.abs(bwidth - kb) / bfactor; //R2B
            } else {
                for (int i = 0; i < r2[0].length; i++) {
                    awidth = peaks.get(2 * i + 1).width50 / (2 * Math.PI);
                    bwidth = peaks.get(2 * i).width50 / (2 * Math.PI);
                    double kex = (awidth + bwidth) / 2;
                    double kb = pb * kex;
                    double ka = (1 - pb) * kex;
                    r2[0][i] = Math.abs(awidth - ka) / afactor; //R2A
                    r2[1][i] = Math.abs(bwidth - kb) / bfactor; //R2B
                }
            }
            return r2;
        } else {
            double[][] r2 = new double[2][1];

            for (int i = 0; i < r2[0].length; i++) {
                double awidth = peaks.get(0).width50 / (2 * Math.PI);
                r2[0][0] = awidth; //R2A
                r2[1][0] = awidth; //R2B
            }
            return r2;
        }
    }

    /**
     * Estimates CEST/R1rho Kex values for initial guesses for fitting.
     *
     * @param peaks List of CESTPeak objects from CEST/R1rho peak guesser.
     * @param fitMode String "cest" or "r1rho" specifying CEST or R1rho data.
     * @return CEST/R1rho Kex value estimate.
     */
    public static double cestKexGuess(List<CESTPeak> peaks, String fitMode) {
        // Estimates CEST kex values from peak widths for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

        double factor = 1;
        if (fitMode.equals("r1rho")) {
            factor = 3;
        }
        if (peaks.size() > 1) {
            double[] kex = new double[peaks.size() / 2];

            for (int i = 0; i < kex.length; i++) {
                double awidth = peaks.get(2 * i + 1).width50 / (2 * Math.PI);
                double bwidth = peaks.get(2 * i).width50 / (2 * Math.PI);
                kex[i] = (awidth + bwidth) / 2; //peaks[2 * i][2]/(2*Math.PI); //Kex
            }
            return kex[0] / factor;
        } else {
            return peaks.get(0).width50 / (2 * Math.PI) / factor; //Kex;
        }

    }

    /**
     * Calculates the CEST/R1rho R1 boundaries for fitting.
     *
     * @param r1 R1 value.
     * @param tex Tex value. Irradiation time.
     * @param delta Delta value that defines the R1 boundary (i.e. R1 +- delta).
     * @return Array with the CEST/R1rho R1 lower and upper bounds.
     */
    public static double[] r1Boundaries(double r1, double tex, double delta) {
        double baseline = Math.exp(-r1 * tex);
        double r1Low = -Math.log(baseline + 0.1) / tex;
        double r1Updiff = baseline - delta;
        if (r1Updiff < 0.01) {
            r1Updiff = 0.01;
        }
        double r1Up = -Math.log(r1Updiff) / tex;
        return new double[]{r1Low, r1Up};
    }

    /**
     * Estimates CEST/R1rho R1 values for initial guesses before fitting.
     *
     * @param yvals Array of the CEST/R1rho intensities.
     * @param Tex Tex value. Irradiation time.
     * @param fitMode String "cest" or "r1rho" specifying CEST or R1rho data.
     * @return Array of CEST/R1rho R1A and R1B value estimates. R1A is listed
     * first, then R1B.
     */
    public static double[] cestR1Guess(double[] yvals, Double Tex, String fitMode) {
        // Estimates CEST R1 values from data baseline intensity and Tex for initial guesses for before fitting.
        // Reference: Palmer, A. G. "Chemical exchange in biomacromolecules: Past, present, and future." J. Mag. Res. 241 (2014) 3-17.

        double[] baseValues = getBaseline(yvals, fitMode);
        double baseline = baseValues[0];
        double[] r1 = {-Math.log(baseline) / 0.3, -Math.log(baseline) / 0.3};
        if (fitMode.equals("r1rho")) {
            r1[0] = baseline; //R1A
            r1[1] = baseline; //R1B
        } else if (fitMode.equals("cest")) {
            if (Tex != null) {
                r1[0] = -Math.log(baseline) / Tex; //R1A
                r1[1] = -Math.log(baseline) / Tex; //R1B
                if (Tex.equals(0.0)) {
                    r1[0] = 0.0;
                    r1[1] = 0.0;
                }
            }
        }
        return r1;
    }

    /**
     * Combines CEST/R1rho x and y values into a single matrix.
     *
     * @param xValues Matrix containing the offset values i.e. CEST/R1rho
     * irradiation frequency (X[0]), the B1 field values (X[1]), and the Tex
     * values (X[2]).
     * @param yValues Array of the CEST/R1rho intensities.
     * @param idNums Array of the ID numbers for the datasets.
     * @param id Integer ID number to retrieve the x and y values.
     * @return Matrix containing the offset (x) and intensity (y) values. X
     * values are in the first array, Y values are in the second.
     */
    public static double[][] getXYValues(double[][] xValues, double[] yValues, int[] idNums, int id) {
        int n = 0;
        for (int idNum : idNums) {
            if (idNum == id) {
                n++;
            }
        }
        double[][] result = new double[xValues.length + 1][n];
        int j = 0;
        for (int i = 0; i < idNums.length; i++) {
            if (idNums[i] == id) {
                for (int k=0;k<xValues.length;k++) {
                    result[k][j] = xValues[k][i];
                }
                result[xValues.length][j] = yValues[i];
                j++;
            }
        }
        return result;

    }

    /**
     * Gets CEST/R1rho indices for a specific dataset.
     *
     * @param idNums Array of the ID numbers for the datasets.
     * @param id Integer ID number to retrieve the x and y values.
     * @return Integer Array containing the dataset indices for the specified ID
     * number.
     */
    public static int[] getIndicies(int[] idNums, int id) {
        int n = 0;
        for (int idNum : idNums) {
            if (idNum == id) {
                n++;
            }
        }
        int[] result = new int[n];
        int j = 0;
        for (int i = 0; i < idNums.length; i++) {
            if (idNums[i] == id) {
                result[j] = i;
                j++;
            }
        }
        return result;

    }

    /**
     * Gets CEST/R1rho x values for a specific dataset.
     *
     * @param xValues Matrix containing the offset values i.e. CEST/R1rho
     * irradiation frequency (X[0]), the B1 field values (X[1]), and the Tex
     * values (X[2]).
     * @param idNums Array of the ID numbers for the datasets.
     * @param id Integer ID number to retrieve the x values.
     * @return Matrix containing the offset (x) values for the specified ID
     * number.
     */
    public static double[][] getXValues(double[][] xValues, int[] idNums, int id) {
        int n = 0;
        for (int idNum : idNums) {
            if (idNum == id) {
                n++;
            }
        }
        double[][] result = new double[xValues.length][n];
        int j = 0;
        for (int i = 0; i < idNums.length; i++) {
            if (idNums[i] == id) {
                for (int k = 0; k < xValues.length; k++) {
                    result[k][j] = xValues[k][i];
                }
                j++;
            }
        }
        return result;

    }

    /**
     * Gets CEST/R1rho y values for a specific dataset.
     *
     * @param values Array of the CEST/R1rho intensities.
     * @param idNums Array of the ID numbers for the datasets.
     * @param id Integer ID number to retrieve the x and y values.
     * @return Matrix containing the intensity (y) values for the specified ID
     * number.
     */
    public static double[] getValues(double[] values, int[] idNums, int id) {
        int n = 0;
        for (int idNum : idNums) {
            if (idNum == id) {
                n++;
            }
        }
        double[] result = new double[n];
        int j = 0;
        for (int i = 0; i < idNums.length; i++) {
            if (idNums[i] == id) {
                result[j] = values[i];
                j++;
            }
        }
        return result;

    }

    /**
     * Purpose - Invoke trained and saved artificial neural network to retrieve
     * guesses for CEST parameters.
     *
     * @param xy - 2D array of x and y values
     * @param peaks - List of peaks returned from cestPeakGuess(...)
     * @param tEx - irradiation time (time of exchange)
     * @param b1Field - radiofrequency field
     *
     * @return guesses - double array containing the neural network's guess for
     * each parameter, in specific order shown below : ORDER = ["Kex", "pB",
     * "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B"]
     * @throws java.io.IOException
     */
    public static double[] cestANNGuess(double[][] xy, List<CESTPeak> peaks, double tEx, double b1Field) throws IOException, Exception {

        double[] xValues = xy[0];
        double[] yValues = xy[1];
        double xMax = xValues[xValues.length - 1];
        double xMin = xValues[0];
        double[] peakScaleVals = {xMax, xMin};
        String[] labelArr = {"kEx", "pB", "deltaA", "deltaB", "r1A", "r1B", "r2A", "r2B"};
        int nPars = labelArr.length; // size : 8
        double[] guesses = new double[nPars];
        HashMap outputs = new HashMap(nPars);

        double baseline = peaks.get(0).baseline;
        // Initialization... Keeping track of guesses
        for (String label : labelArr) {
            outputs.put(label, 0.0);
        }

        // Loading ANN ...
        ANNLoader ANN = ANNLoader.getInstance("data/ANNR1RhoPerturbation_A.txt"); // for B, depth is scaled b/t 1,-1
        ArtificialNeuralNetwork trainedNetwork = ANN.getTrainedNetwork();
        if (trainedNetwork != null) {
            // updating DELTA_A and DELTA_B
            outputs.replace("deltaA", peaks.get(0).getPosition());
            outputs.replace("deltaB", peaks.get(peaks.size() - 1).getPosition());

            // Get scale values
            HashMap scaleTracker = ANN.getScaleValues();

            // updating scale values for deltaA and deltaB
            scaleTracker.replace("deltaA", peakScaleVals);
            scaleTracker.replace("deltaB", peakScaleVals);

            // CODE TO SCALE INPUT:
            List<Double> annInput = new ArrayList(ANN.getNumberOfInputNodes());
            peaks.forEach((peak) -> {
                annInput.add(peak.getDepth());
                for (double width : peak.getWidths()) {
                    annInput.add(Utilities.scale(width, (double[]) scaleTracker.get("b1Field"), true));
                }
            });

            // Guessing R1A/R1B
            double[] r1 = cestR1Guess(yValues, tEx, "cest");
            outputs.replace("r1A", r1[0]);
            outputs.replace("r1B", r1[1]);
            annInput.add(Utilities.scale(b1Field, (double[]) scaleTracker.get("b1Field"), true));
            annInput.add(Utilities.scale(tEx, (double[]) scaleTracker.get("Tex"), true));
            annInput.add(baseline);

            // CODE TO INVOKE NEURAL NETWORK:
            Access1D wrappedInput = Access1D.wrap(annInput);
            MatrixStore<Double> result = trainedNetwork.invoke(wrappedInput);
            double[] resultArr = result.toRawCopy1D();

            // CODE TO CONVERT NEURAL NETWORK OUTPUT:
            String[] networkLabelArr = {"kEx", "pB", "r2A", "r2B"};
            double[] unscaledResultArr = new double[resultArr.length];

            for (int l = 0; l < resultArr.length; l++) {
                if (scaleTracker.containsKey(networkLabelArr[l])) {
                    unscaledResultArr[l] = Utilities.scale(resultArr[l], (double[]) scaleTracker.get(networkLabelArr[l]), false);
                    outputs.replace(networkLabelArr[l], unscaledResultArr[l]);
                }
            }
            for (int j = 0; j < labelArr.length; j++) {
                guesses[j] = (double) outputs.get(labelArr[j]);
            }
            return guesses; // should be 8, ordered properly
        } else {
            throw new Exception("Could not load neural network.");
        }
    }

}
