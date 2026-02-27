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
package org.comdnmr.modelfree;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static org.comdnmr.modelfree.RelaxFit.DiffusionType;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.ISOTROPIC;

/**
 *
 * @author brucejohnson
 */
public class RelaxEquations {

    static Map<String, RelaxEquations> relaxMap = new HashMap<>();

    static final int S = 1;
    static final int ImS = 2;
    static final int I = 3;
    static final int IpS = 4;
    public final static double SQRT2 = Math.sqrt(2.0);

    public final static double MU0 = 4.0e-7 * Math.PI;
    public final static double GAMMA_N = -2.7116e7;
    public final static double GAMMA_C = 6.72828e7;
    public final static double GAMMA_H = 2.6752218744e8;
    public final static double GAMMA_D = 4.1065e7;

    public final static double HBAR = 1.054571817e-34;
    public final static double R_HH = 3.0e-10; // Simon chucked this in in order for some tests to work.
    public final static double R_HN = 1.02e-10;
    public final static double R_HC = 1.09e-10;
    public final static double R_CC = 1.51e-10;
    public final static double SIGMA = -172.0e-6;  // N.B. this is for Δσ = (3/2)σ = σzz - (1/2) (σxx + σyy)
    public final static double QCC = Math.PI*167.0e3/2.0;
    public final static double QCC2 = QCC * QCC;

    public final static Map<String, Double> GAMMA_MAP = new HashMap<>();
    public final static Map<String, Double> R_MAP = new HashMap<>();
    public final static Map<String, Double> SIGMA_MAP = new HashMap<>();

    static {
        GAMMA_MAP.put("H", GAMMA_H);
        GAMMA_MAP.put("N", GAMMA_N);
        GAMMA_MAP.put("1H", GAMMA_H);
        GAMMA_MAP.put("15N", GAMMA_N);
        GAMMA_MAP.put("C", GAMMA_C);
        GAMMA_MAP.put("13C", GAMMA_C);
        GAMMA_MAP.put("D", GAMMA_D);
        R_MAP.put("HH", R_HH);  // Simon chucked this in in order for some tests to work.
        R_MAP.put("HN", R_HN);
        R_MAP.put("NH", R_HN);
        R_MAP.put("HC", R_HC);
        R_MAP.put("CH", R_HC);
        R_MAP.put("DC", R_HC);
        R_MAP.put("CD", R_HC);
        SIGMA_MAP.put("H", 0.0);  // Simon put this in to support R1RhoAA
        SIGMA_MAP.put("N", SIGMA);
        SIGMA_MAP.put("C", SIGMA);
    }

    private final double r;
    private final double d;
    private final double d2;
    private final double c;  // (ωΔσ) / √3
    private final double c2;   // (ω²Δσ²) / 3

    private final double gammaS;
    private final double gammaI;
    private final double sf;
    private final double wI;
    private final double wS;
    double[] wValues;

    //   consider using scaled versions (smaller exponents)
    /**
     *
     * @param sf double. 1H NMR Spectrometer frequency (MHz).
     * @param elem1 String. First element ("H" for 1H NMR).
     * @param elem2 String. Second element ("C", "N", etc.)
     */
    public RelaxEquations(double sf, String elem1, String elem2) {
        gammaI = GAMMA_MAP.get(elem1);
        gammaS = GAMMA_MAP.get(elem2);
        if (elem1.equals("D")) {
            wI = sf * 2.0 * Math.PI * gammaI / GAMMA_H;
            wS = sf * 2.0 * Math.PI * gammaS / GAMMA_H;
        } else {
            wI = sf * 2.0 * Math.PI;
            wS = wI * gammaS / gammaI;
        }
        r = R_MAP.get(elem1 + elem2);
        d = -MU0 * (gammaI * gammaS * HBAR) / (4.0 * Math.PI * Math.pow(r, 3.0));
        d2 = d * d;
        c = wS * SIGMA_MAP.get(elem2) / Math.sqrt(3.0);
        c2 = c * c;
        if (elem1.equals("D")) {
            wValues = new double[]{0.0, wI, 2.0 * wI};
        } else {
            wValues = new double[]{0.0, wS, wI - wS, wI, wI + wS};
        }
        this.sf = sf;
    }

    public static RelaxEquations getRelaxEquations(double sf, String elem1, String elem2) {
        int sfI = (int) Math.round(sf / 1.0e6);
        String key = sfI + elem1 + elem2;
        if (!relaxMap.containsKey(key)) {
            RelaxEquations rObj = new RelaxEquations(sf, elem1, elem2);
            relaxMap.put(key, rObj);
        }
        return relaxMap.get(key);
    }

    public static void setR(String elem1, String elem2, double value) {
        R_MAP.put(elem1+elem2, value);
        R_MAP.put(elem2+elem1, value);
        relaxMap.clear();
    }
    public static double getR(String elem1, String elem2) {
        return R_MAP.get(elem1+elem2);
    }

    public static void setSigma(String elem1, double value) {
        SIGMA_MAP.put(elem1, value);
        relaxMap.clear();
    }
    public static double getSigma(String elem1) {
        return SIGMA_MAP.get(elem1);
    }

    public double getSF() {
        return sf;
    }

    public double getGammaS() {
        return gammaS;
    }

    public double getGammaI() {
        return gammaI;
    }

    public double getD() {
        return d;
    }

    public double getD2() {
        return d2;
    }

    public double getC2() {
        return c2;
    }

    public double getR() {
        return r;
    }

    public double[] getW() {
        return wValues;
    }

    public double getWI() {
        return wI;
    }

    public double getWS() {
        return wS;
    }

    public static double getSF(double sf, String elemX) {
        return Math.abs(sf * GAMMA_MAP.get(elemX) / GAMMA_H);
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax.
    /**
     * Model Free spectral density function, J(omega), calculation using Model
     * 1.
     *
     * @param w double. The frequency, omega.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @return J(w) value.
     */
    public double JModelFree(double w, double tauM, double s2) {
        double value1 = s2 / (1.0 + w * w * tauM * tauM);
        return 0.4 * tauM * (value1);
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax.
    // tau = ts in Art Palmer's code (taue in the paper: Phys Chem Chem Phys, 2016, 18, 5839-5849), and taue in Relax.
    /**
     * Model Free spectral density function, J(omega), calculation using Model
     * 2.
     *
     * @param w double. The frequency, omega.
     * @param tau double. The internal correlation time.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @return J(w) value.
     */
    public double JModelFree(double w, double tau, double tauM, double s2) {
        double value1 = s2 / (1.0 + w * w * tauM * tauM);
        double value2 = ((1.0 - s2) * (tau + tauM) * tau) / ((tau + tauM) * (tau + tauM) + w * w * tauM * tauM * tau * tau);
        return 0.4 * tauM * (value1 + value2);
    }

    // Note: tauM = tm in Art Palmer's code. tau = ts in Art Palmer's code.
    /**
     * Model Free spectral density function, J(omega), calculation using Model
     * 5.
     *
     * @param w double. The frequency, omega.
     * @param tau double. The internal correlation time.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @param sf2 double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @return J(w) value.
     */
    public double JModelFree(double w, double tau, double tauM, double s2, double sf2) {
        double value1 = s2 / (1.0 + w * w * tauM * tauM);
        double value2 = ((sf2 - s2) * (tau + tauM) * tau) / ((tau + tauM) * (tau + tauM) + w * w * tauM * tauM * tau * tau);
        return 0.4 * tauM * (value1 + value2);
    }

    /**
     * Model Free spectral density function, J(omega), calculation using Model
     * 6.
     *
     * @param w double. The frequency, omega.
     * @param tauF double. The internal fast correlation time.
     * @param tauM double. The overall correlation time.
     * @param tauS double. The internal slow correlation time.
     * @param s2 double. The order parameter S^2.
     * @param sf2 double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @return J(w) value.
     */
    public double JModelFree(double w, double tauF, double tauM, double tauS, double s2, double sf2) {
        double value1 = s2 / (1.0 + w * w * tauM * tauM);
        double value2 = ((1 - sf2) * (tauF + tauM) * tauF) / ((tauF + tauM) * (tauF + tauM) + w * w * tauM * tauM * tauF * tauF);
        double value3 = ((sf2 - s2) * (tauS + tauM) * tauS) / ((tauS + tauM) * (tauS + tauM) + w * w * tauM * tauM * tauS * tauS);
        return 0.4 * tauM * (value1 + value2 + value3);
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax.
    /**
     * Model Free spectral density function, J(omega), calculations using Model
     * 1.
     *
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJModelFree(double tauM, double s2) {
        double J0 = JModelFree(0.0, tauM, s2);
        double JIplusS = JModelFree(wI + wS, tauM, s2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tauM, s2); // R1, R2, NOE
        double JI = JModelFree(wI, tauM, s2); // R2
        double JS = JModelFree(wS, tauM, s2); //R1, R2
        return new double[]{J0, JS, JIminusS, JI, JIplusS};
    }

    public double[] getOmegas() {
        double[] result = new double[5];
        result[0] = 0.0;
        result[S] = wS;
        result[ImS] = wI-wS;
        result[I] = wI;
        result[IpS] = wI + wS;
        return result;
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax. tau = ts in Art Palmer's code (taue in the paper), and taue in Relax.
    /**
     * Model Free spectral density function, J(omega), calculations using Model
     * 2.
     *
     * @param tau double. The internal correlation time.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJModelFree(double tau, double tauM, double s2) {
        double J0 = JModelFree(0.0, tau, tauM, s2);
        double JIplusS = JModelFree(wI + wS, tau, tauM, s2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tau, tauM, s2); // R1, R2, NOE
        double JI = JModelFree(wI, tau, tauM, s2); // R2
        double JS = JModelFree(wS, tau, tauM, s2); //R1, R2
        return new double[]{J0, JS, JIminusS, JI, JIplusS};
    }

    // Note: tauM = tm in Art Palmer's code. tau = ts in Art Palmer's code.
    /**
     * Model Free spectral density function, J(omega), calculations using Model
     * 5.
     *
     * @param tau double. The internal correlation time.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @param sf2 double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJModelFree(double tau, double tauM, double s2, double sf2) {
        double J0 = JModelFree(0.0, tau, tauM, s2, sf2);
        double JIplusS = JModelFree(wI + wS, tau, tauM, s2, sf2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tau, tauM, s2, sf2); // R1, R2, NOE
        double JI = JModelFree(wI, tau, tauM, s2, sf2); // R2
        double JS = JModelFree(wS, tau, tauM, s2, sf2); //R1, R2
        return new double[]{J0, JS, JIminusS, JI, JIplusS};
    }

    /**
     * Model Free spectral density function, J(omega), calculations using Model
     * 6.
     *
     * @param tauF double. The internal fast correlation time.
     * @param tauM double. The overall correlation time.
     * @param tauS double. The internal slow correlation time.
     * @param s2 double. The order parameter S^2.
     * @param sf2 double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJModelFree(double tauF, double tauM, double tauS, double s2, double sf2) {
        double J0 = JModelFree(0.0, tauF, tauM, tauS, s2, sf2);
        double JIplusS = JModelFree(wI + wS, tauF, tauM, tauS, s2, sf2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tauF, tauM, tauS, s2, sf2); // R1, R2, NOE
        double JI = JModelFree(wI, tauF, tauM, tauS, s2, sf2); // R2
        double JS = JModelFree(wS, tauF, tauM, tauS, s2, sf2); //R1, R2
        return new double[]{J0, JS, JIminusS, JI, JIplusS};
    }

    /**
     * Model Free spectral density function, J(omega), diffusion calculation
     * using ModelFree Model 1, 2, 5, or 6. Equations from the SI of Berlin K.;
     * Longhini, A.; Dayie, T. K. and Fushman, D., J. Biomol NMR, 2013.
     *
     * @param diffType DiffusionType. The type of diffusion: anisotropic,
     * prolate, or oblate.
     * @param w double. The frequency, omega.
     * @param D double[][]. The diagonalized diffusion matrix.
     * @param VT double[][]. The transposed orthonormal matrix of the
     * eigenvectors of D.
     * @param v double[]. The unit vector for the SI bond.
     * @param s2 double. The order parameter S^2.
     * @param tauF Double. The internal fast correlation time.
     * @param sf2 Double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @param tauS Double. The internal slow correlation time.
     * @return J(w) value.
     */
    private double JDiffusion(DiffusionType diffType, double w, double[][] D, double[][] VT, double[] v, double s2, Double tauF, Double sf2, Double tauS) {
        double w2 = w * w;
        double sum = 0.0;
        if (diffType == ISOTROPIC) {
            double Diso = (D[0][0] + D[1][1] + D[2][2]) / 3.0;
            double tauC = 1.0 / (6 * Diso);
            double value1 = s2 * tauC / (1.0 + w2 * tauC * tauC);
            double value2 = 0.0;
            double value3 = 0.0;
            if (tauF != null && sf2 == null && tauS == null) {
                double tauf = tauC * tauF / (tauC + tauF);
                value2 = (1.0 - s2) * (tauf) / (1.0 + w2 * tauf * tauf);
            } else if (tauF != null && sf2 != null && tauS == null) {
                double tauf = tauC * tauF / (tauC + tauF);
                value2 = (sf2 - s2) * (tauf) / (1.0 + w2 * tauf * tauf);
            } else if (tauF != null && sf2 != null && tauS != null) {
                double taue = tauC * tauS / (tauC + tauS);
                double tauf = tauC * tauF / (tauC + tauF);
                value2 = (sf2 - s2) * (taue) / (1.0 + w2 * taue * taue);//(eS[i]*eS[i] + w2);
                value3 = (1.0 - sf2) * (tauf) / (1.0 + w2 * tauf * tauf);//(eF[i]*eF[i] + w2);
            }
            sum = value1 + value2 + value3;
        } else {
            DiffusionPars diffPars = new DiffusionPars(diffType, D, VT, v);
            double[] dDiff = diffPars.dDiff;
            double[] a = diffPars.a;
            double[] eF = null;
            double[] eS = null;
            if (tauF != null) {
                eF = diffPars.calcDiffusione(tauF);
                if (tauS != null) {
                    eS = diffPars.calcDiffusione(tauS);
                }
            }
            double[] Df = diffPars.getDf(w2);
            for (int i = 0; i < dDiff.length; i++) {
                double value1 = s2 * (Df[i] * a[i]);
                double value2 = 0.0;
                double value3 = 0.0;
                if (tauF != null && sf2 == null && tauS == null) {
                    value2 = (1.0 - s2) * (eF[i] * a[i]) / (1.0 + w2 * eF[i] * eF[i]);//(eF[i]*eF[i] + w2);
                } else if (tauF != null && sf2 != null && tauS == null) {
                    value2 = (sf2 - s2) * (eF[i] * a[i]) / (1.0 + w2 * eF[i] * eF[i]);//(eF[i]*eF[i] + w2);
                } else if (tauF != null && sf2 != null && tauS != null) {
                    value2 = (sf2 - s2) * (eS[i] * a[i]) / (1.0 + w2 * eS[i] * eS[i]);//(eS[i]*eS[i] + w2);
                    value3 = (1.0 - sf2) * (eF[i] * a[i]) / (1.0 + w2 * eF[i] * eF[i]);//(eF[i]*eF[i] + w2);
                }
                sum += value1 + value2 + value3;
            }
        }
        return 0.4 * sum;
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax.
    /**
     * Model Free spectral density function, J(omega), calculations using
     * ModelFree Model 1, 2, 5, or 6.
     *
     * @param diffType DiffusionType. The type of diffusion: anisotropic,
     * prolate, or oblate.
     * @param D double[][]. The diagonalized diffusion matrix.
     * @param VT double[][]. The transposed orthonormal matrix of the
     * eigenvectors of D.
     * @param v double[]. The unit vector for the SI bond.
     * @param s2 double. The order parameter S^2.
     * @param tauF Double. The internal fast correlation time.
     * @param sf2 Double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @param tauS Double. The internal slow correlation time.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJDiffusion(DiffusionType diffType, double[][] D, double[][] VT, double[] v, double s2, Double tauF, Double sf2, Double tauS) {
        double J0 = JDiffusion(diffType, 0.0, D, VT, v, s2, tauF, sf2, tauS);
        double JIplusS = JDiffusion(diffType, wI + wS, D, VT, v, s2, tauF, sf2, tauS); //R1, R2, NOE
        double JIminusS = JDiffusion(diffType, wI - wS, D, VT, v, s2, tauF, sf2, tauS); // R1, R2, NOE
        double JI = JDiffusion(diffType, wI, D, VT, v, s2, tauF, sf2, tauS); // R2
        double JS = JDiffusion(diffType, wS, D, VT, v, s2, tauF, sf2, tauS); //R1, R2
        return new double[]{J0, JS, JIminusS, JI, JIplusS};
    }

    /**
     * Spectral density function calculation.
     *
     * @param w double. The frequency, omega.
     * @param tau double. The correlation time.
     * @return double. The spectral density, J(w).
     */
    public double J(double w, double tau) {
        return 0.4 * tau / (1.0 + w * w * tau * tau);
    }

    /**
     * Spectral density function calculation.
     *
     * @param w double. The frequency, omega.
     * @param tau double. The correlation time.
     * @param S2 double. The order parameter.
     * @return double. The spectral density, J(w).
     */
    public double J(double w, double tau, double S2) {
        return 0.4 * S2 * tau / (1.0 + w * w * tau * tau);
    }

    /*
     * Model-free function for CH3 groups. tf is time scale for methyl rotation
     * Derived from Skrynnikov and Kay for Sf2 = 1
     *
     */
    public double getJCH3(double w, double tauM, double Sf2, double Ss2, double tauF, double tauS) {
        double Sf2p = Sf2 / 9.0;
        double t1 = tauM * tauS / (tauM + tauS);
        double t2 = tauM * tauF / (tauM + tauF);
        double t3 = tauM * tauS * tauF / (tauM * tauS + tauM * tauF + tauS * tauF);
        return (2.0 / 5.0) * (Sf2p * Ss2 * tauM / (1 + w * w * tauM * tauM) +
                Sf2p * (1 - Ss2) * t1 / (1 + w * w * t1 * t1) +
                (1 - Sf2p) * Ss2 * t2 / (1 + w * w * t2 * t2) +
                (1 - Sf2p) * (1 - Ss2) * t3 / (1 + w * w * t3 * t3));
    }

    public double getJCH3red(double w, double tauM, double Sf2, double tauF) {
        double Sf2p = Sf2 / 9.0;
        double t1 = tauM * tauF / (tauM + tauF);
        return (2.0 / 5.0) * (Sf2p * tauM / (1 + w * w * tauM * tauM) +
                (1 - Sf2p) * t1 / (1 + w * w * t1 * t1));
    }

    public double getJCH3min(double w, double tauM, double Sf2) {
        double Sf2p = Sf2 / 9.0;
        return (2.0 / 5.0) * (Sf2p * tauM / (1.0 + w * w * tauM * tauM));
    }

    public double[] getDiffusionConstants(String type) {
        String[] types = {"sphere", "spheroid", "ellipsoid"};
        double[][] constants = {{1}};//,
//            {0.25*(3.0*dz2 - 1)*(3.0*dz2 - 1), 3*dz2*(1 - dz2), 0.75*(3.0*dz2 - 1)*(3.0*dz2 - 1)},
//            {0.25*(dtot - e), 3*dy2*dz2, 3*dx2*dz2, 3*dx2*dy2, 0.25*(dtot + e)}};
        return constants[Arrays.asList(types).indexOf(type)];
    }

    public double[] getCorrelationTimes(String type, double Diso, double Da, double Dr) {
        String[] types = {"sphere", "spheroid", "ellipsoid"};
        double[][] tauInv = {{6 * Diso}};//,
//            {6*Diso - 2*Da, 6*Diso - Da, 6*Diso + 2*Da},
//            {6*Diso - 2*Da*R, 6*Diso - Da*(1 + 3*Dr), 6*Diso - Da*(1 - 3*Dr), 6*Diso + 2*Da, 6*Diso + 2*Da*R}};
        int index = Arrays.asList(types).indexOf(type);
        double[] taus = new double[tauInv[index].length];
        for (int i = 0; i < taus.length; i++) {
            taus[i] = 1 / tauInv[index][i];
        }
        return taus;
    }

    public double calcS2(double J0, double bN, double mN) {
        return (5.0 / 2) * Math.sqrt((J0 - bN) * mN);
    }

    /**
     * Spectral density function calculations.
     *
     * @param tau double. The correlation time.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJ(double tau) {
        double J0 = J(0.0, tau); // R2
        double JIplusS = J(wI + wS, tau); //R1, R2, NOE
        double JIminusS = J(wI - wS, tau); // R1, R2, NOE
        double JI = J(wI, tau); // R2
        double JS = J(wS, tau); //R1, R2
        return new double[]{J0, JS, JIminusS, JI, JIplusS};
    }

    /**
     * Spectral density function calculations.
     *
     * @param tau double. The correlation time.
     * @param S double. The order parameter.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJ(double tau, double S) {
        double[] result = getJ(tau);
        for (int i = 0; i < result.length; i++) {
            result[i] *= S;
        }
        return result;
    }

    /**
     * R2:R1 ratio calculation.
     *
     * @param tau double. The correlation time.
     * @return double. R2/R1 value.
     */
    public double r2r1Ratio(double tau) {
        double num = 4.0 * J(0, tau) + J(wS - wI, tau) + 3.0 * J(wS, tau)
                + 6.0 * J(wI, tau) + 6.0 * J(wS + wI, tau)
                + (c2 / (3.0 * d2) * (4.0 * J(0, tau) + 3.0 * J(wS, tau)));

        double denom = 2.0 * J(wS - wI, tau) + 6.0 * J(wS, tau) + 12.0 * J(wS + wI, tau) + 2.0 * (c2 / (3.0 * d2) * J(wS, tau));
        return num / denom;
    }

    /**
     * R1 calculation
     *
     * @param tau double. The correlation time.
     * @return double. R1 value.
     */
    public double R1(double tau) {
        double dipolarContrib = d2 / 4.0 * (J(wI - wS, tau) + 3.0 * J(wS, tau) + 6.0 * J(wI + wS, tau));
        double csaContrib = c2 * J(wS, tau);
        return dipolarContrib + csaContrib;
    }

    /**
     * R2 calculation
     *
     * @param tau double. The correlation time.
     * @return double. R2 value.
     */
    public double R2(double tau) {
        double dipolarContrib = d2 / 8.0 * (4.0 * J(0.0, tau) + J(wI - wS, tau) + 3.0 * J(wS, tau)
                + 6.0 * J(wI, tau) + 6.0 * J(wI + wS, tau));
        double csaContrib = c2 / 6 * (4.0 * J(0.0, tau) + 3.0 * J(wS, tau));
        return dipolarContrib + csaContrib;
    }

    /**
     * NOE calculation
     *
     * @param tau double. The correlation time.
     * @return double. NOE value.
     */
    public double NOE(double tau) {
        double R1 = R1(tau);
        return 1.0 + (d2 / (4.0 * R1)) * (gammaI / gammaS)
                * (6.0 * J(wI + wS, tau) - J(wI - wS, tau));
    }

    /**
     * R1 calculation
     *
     * @param J double[]. Array of spectral density function values, J(w).
     * @return double. R1 value.
     */
    public double R1(double[] J) {
        double dipolarContrib = d2 / 4.0 * (J[ImS] + 3.0 * J[S]
                + 6.0 * J[IpS]);
        double csaContrib = c2 * J[S];
        return dipolarContrib + csaContrib;
    }

    /**
     * R2 calculation
     *
     * @param J double[]. Array of spectral density function values, J(w).
     * @param Rex double. Rate of exchange, Rex, value.
     * @return double. R2 value.
     */
    public double R2(double[] J, double Rex) {
        double dipolarContrib = d2 / 8.0 * (4.0 * J[0] + J[ImS]
                + 3.0 * J[S]
                + 6.0 * J[I] + 6.0 * J[IpS]);
        double csaContrib = c2 / 6 * (4.0 * J[0] + 3.0 * J[S]);
        return dipolarContrib + csaContrib + Rex;
    }

    public double R1_D(double[] J) {
        double jw = J[1];
        double j2w = J[2];
        return 3.0 * QCC2 * (jw + 4.0 * j2w);
    }

    public double R2_D(double[] J) {
        double j0 = J[0];
        double jw = J[1];
        double j2w = J[2];
        return (3.0 / 2.0) * QCC2 * (3.0 * j0 + 5.0 * jw + 2.0 * j2w);
    }

    public double RQ_D(double[] J) {
        double jw = J[1];
        return 9.0 * QCC2 * jw;
    }

    public double Rap_D(double[] J) {
        double j0 = J[0];
        double jw = J[1];
        double j2w = J[2];
        return (3.0 / 2.0) * QCC2 * (3 * j0 + jw + 2 * j2w);
    }

    /**
     * NOE calculation
     *
     * @param J double[]. Array of spectral density function values.
     * @return double. NOE value.
     */
    public double NOE(double[] J) {
        double R1 = R1(J);
        return 1.0 + (d2 / (4.0 * R1)) * (gammaI / gammaS)
                * (6.0 * J[IpS] - J[ImS]);
    }

    /**
     * Sigma calculation
     *
     * @param J double[]. Array of spectral density function values.
     * @return double. Sigma value.
     */
    public double sigmaSI(double[] J) {
        double R1 = R1(J);
        double NOE = NOE(J);
        return R1 * (NOE - 1) * (gammaS / gammaI);
    }

    /**
     * Gamma calculation
     *
     * @param J double[]. Array of spectral density function values, J(w).
     * @param Rex double. Rate of exchange, Rex, value.
     * @return double. Gamma value.
     */
    public double gamma(double[] J, double Rex) {
        double R1 = R1(J);
        double R2 = R2(J, Rex);
        double sigmaSI = sigmaSI(J);
        return R2 - 0.5 * R1 - 0.454 * sigmaSI;
    }

    public double TRACTdeltaAlphaBeta(double tauC) {
        double B0 = wI / GAMMA_H;

        double ddN = 160.0e-6;
        double theta = 17.0 * Math.PI / 180.0;
        double p = MU0 * gammaI * gammaS * HBAR / (8.0 * Math.PI * SQRT2 * r * r * r);
        double dN = gammaS * B0 * ddN / (3.0 * SQRT2);

        double cosTheta = Math.cos(theta);
        double[] J = getJ(tauC);
        return 2.0 * p * dN * (4.0 * J[0] + 3.0 * J[S]) * (3.0 * cosTheta * cosTheta - 1.0);
    }

    /**
     * Calculate rhoExp from Eqn 6 in Berlin K.; Longhini, A.; Dayie, T. K. and
     * Fushman, D., J. Biomol NMR, 2013.
     *
     * @param R1 double. Experimentally measured R1 value.
     * @param R2 double. Experimentally measured R2 value.
     * @param NOE double. Experimentally measured NOE value.
     * @param J double[]. Array of spectral density function values, J(w).
     * @return double. RhoExp.
     */
    public double calcRhoExp(double R1, double R2, double NOE, double[] J) {
        double w1 = (J[ImS] + 6 * J[IpS]) / (6 * J[IpS] - J[ImS]);
        double w2 = (6 * J[I]) / (6 * J[IpS] - J[ImS]);
        double f = (gammaS / gammaI) * R1 * (NOE - 1);
        return (2 * R2 - R1 - w2 * f) / (R1 - w1 * f);
    }

    /**
     * Calculate the standard deviation error in rhoExp.
     *
     * @param R1 double. Experimentally measured R1 value.
     * @param R2 double. Experimentally measured R2 value.
     * @param NOE double. Experimentally measured NOE value.
     * @param J double[]. Array of spectral density function values, J(w).
     * @param R1err double. Error value for the experimentally measured R1
     * value.
     * @param R2err double. Error value for the experimentally measured R2
     * value.
     * @param NOEerr double. Error value for the experimentally measured NOE
     * value.
     * @param rhoExp double. Calculated rhoExp value.
     * @return double. RhoExp sigma squared error.
     */
    public double calcRhoExpError(double R1, double R2, double NOE, double[] J, double R1err, double R2err, double NOEerr, double rhoExp) {
        double w1 = (J[ImS] + 6 * J[IpS]) / (6 * J[IpS] - J[ImS]);
        double w2 = (6 * J[I]) / (6 * J[IpS] - J[ImS]);
        double f = (gammaS / gammaI) * R1 * (NOE - 1);
        double rhoExp1 = 2 * R2 - R1 - w2 * f;
        double rhoExp2 = R1 - w1 * f;
        double v = (R1err / R1) * (R1err / R1) + (NOEerr / NOE) * (NOEerr / NOE);
        double rhoExpErr1 = Math.sqrt(2 * 2 * R2err * R2err + R1err * R1err
                + w2 * w2 * (f * f * v + R1err * R1err));
        double rhoExpErr2 = Math.sqrt(R1err * R1err
                + w1 * w1 * (f * f * v + R1err * R1err));
        return Math.sqrt(rhoExp * rhoExp * ((rhoExpErr1 / rhoExp1) * (rhoExpErr1 / rhoExp1) + (rhoExpErr2 / rhoExp2) * (rhoExpErr2 / rhoExp2)));
    }

    /**
     * Calculate rhoExp from function computeRho() in RotDif code
     * (..../relax/RelaxationDatum.java).
     *
     * @param R1 double. Experimentally measured R1 value.
     * @param R2 double. Experimentally measured R2 value.
     * @param NOE double. Experimentally measured NOE value.
     * @param J double[]. Array of spectral density function values, J(w).
     * @return double. RhoExp.
     */
    public double calcRhoExpCode(double R1, double R2, double NOE, double[] J) {
        double w1 = (J[ImS] + 6 * J[IpS]) / (6 * J[IpS] - J[ImS]);
        double w2 = 0.5 * (6 * J[I] + 6 * J[IpS] + J[ImS]) / (6 * J[IpS] - J[ImS]);
        double f = (gammaS / gammaI) * R1 * (NOE - 1);
        double R1c = R1 - w1 * f;
        double R2c = R2 - w2 * f;
        return (2 * R2c - R1c) / (R1c);
    }

    /**
     * Calculate rhoPred from Eqn 7 in Berlin K.; Longhini, A.; Dayie, T. K. and
     * Fushman, D., J. Biomol NMR, 2013.
     *
     * @param J double[]. Array of spectral density function values, J(w).
     * @return double. RhoExp.
     */
    public double calcRhoPred(double[] J) {
        return (4.0 / 3.0) * (J[0] / J[S]);
    }

    // SOLID STATE STUFF

    /**
    * Compute the CSA contribution to R1rho relaxation in the solid-state context.
    * <p>
    * Computes Eq. 39 in J. Phys. Chem. B 2017, 121, 25, 6117–6130.
    *
    * @param wr Rotating frame (MAS) frequency (rad s-1).
    * @param we Spin-lock field frequency (rad s-1).
    * @param tau Rotational correlation time (s).
    * @param S2 Order parameter.
    * @return R1rhoCSA: Contribution from CSA to R1rho.
    */
    public double r1RhoCSA(double wr, double we, double tau, double S2) {
        // Using public double J(double w, double tau, double S2)
        double oneMinusS2 = 1.0 - S2;
        // Usually, R1rho experiments are performed with the locking field on resonance,
        // such that theta_p is 90.
        // In general, arccos(2 * pi * offset / we) should be used.
        double theta_p = 0.5 * Math.PI;
        // N.B. c2 is equivalent to (ω² Δσ²) / 3
        return c2 * (
            (1.0 / 9.0) * Math.pow(Math.sin(theta_p), 2.0) * (
                J(we - 2 * wr, tau, oneMinusS2) +
                2.0 * J(we - wr, tau, oneMinusS2) +
                2.0 * J(we + wr, tau, oneMinusS2) +
                J(we + 2 * wr, tau, oneMinusS2)
            ) +
            0.25 * (3.0 + Math.cos(2 * theta_p)) * J(wS, tau, oneMinusS2)
        );
    }

    /**
     * Compute the heteronuclear dipolar contribution to R1rho relaxation in the
     * solid-state context.
     * <p>
     * Computes Eq. 33 in J. Phys. Chem. B 2017, 121, 25, 6117–6130.
     *
     * @param wr Rotating frame (MAS) frequency (rad s-1).
     * @param we Spin-lock field frequency (rad s-1).
     * @param tau Rotational correlation time (s).
     * @param S2 Order parameter.
     * @return R1rhoDD: Contribution from dipolar interaction to R1rho.
     */
    public double r1RhoIS(double wr, double we, double tau, double S2) {
        // See applicable comments on omeMinusS2 and theta_p in R1rhoCSA
        double oneMinusS2 = 1.0 - S2;
        double theta_p = 0.5 * Math.PI;
        return 0.25 * d2 * (
            (1.0 / 3.0) * Math.pow(Math.sin(theta_p), 2.0) * (
                J(we - 2 * wr, tau, oneMinusS2) +
                2.0 * J(we - wr, tau, oneMinusS2) +
                2.0 * J(we + wr, tau, oneMinusS2) +
                J(we + 2 * wr, tau, oneMinusS2) +
                9.0 * J(wS, tau, oneMinusS2)
            ) +
            0.25 * (3.0 + Math.cos(2 * theta_p)) * (
                3.0 * J(wI, tau, oneMinusS2) +
                6.0 * J(wI + wS, tau, oneMinusS2) +
                J(wI - wS, tau, oneMinusS2)
            )
        );
    }

    public double r1RhoAA(double wr, double we, double tau, double S2) {
        if (getGammaI() != getGammaS()) {
            throw new RuntimeException("Expected the two nuclei to be the same.");
        }
        // See applicable comments on omeMinusS2 and theta_p in R1rhoCSA
        double oneMinusS2 = 1.0 - S2;
        double theta_p = 0.5 * Math.PI;
        double sinSq2theta = Math.pow(Math.sin(2.0 * theta_p), 2.0);
        double sinFourthTheta = Math.pow(Math.sin(theta_p), 4.0);
        double threeCos2theta = 3.0 * Math.cos(2.0 * theta_p);
        double r1Rho = 0.75 * d2 * (
            0.125 * sinSq2theta * (
                J(wr - we, tau, oneMinusS2) +
                J(wr + we, tau, oneMinusS2)
            ) +
            0.5 * sinFourthTheta * (
                J(wr - 2.0 * we, tau, oneMinusS2) +
                J(wr + 2.0 * we, tau, oneMinusS2)
            ) +
            0.0625 * sinSq2theta * (
                J(2.0 * wr - we, tau, oneMinusS2) +
                J(2.0 * wr + we, tau, oneMinusS2)
            ) +
            0.25 * sinFourthTheta * (
                J(2.0 * wr - 2.0 * we, tau, oneMinusS2) +
                J(2.0 * wr + 2.0 * we, tau, oneMinusS2)
            ) +
            0.25 * (7.0 - threeCos2theta) * J(wI, tau, oneMinusS2) +
            0.5 * (5.0 + threeCos2theta) * J(2.0 * wI, tau, oneMinusS2)
        );
        return r1Rho;
    }

    public double r1RhoAB(double wr, double we, double tau, double S2) {
        if (getGammaI() != getGammaS()) {
            throw new RuntimeException("Expected the two nuclei to be the same.");
        }
        // See applicable comments on omeMinusS2 and theta_p in R1rhoCSA
        double oneMinusS2 = 1.0 - S2;
        double theta_p = 0.5 * Math.PI;
        double cos2theta = Math.cos(2.0 * theta_p);
        double cosSq2theta = Math.pow(cos2theta, 2.0);
        double sinSq2theta = Math.pow(Math.sin(2.0 * theta_p), 2.0);
        double sinFourthTheta = Math.pow(Math.sin(theta_p), 4.0);
        double r1Rho = 0.25 * d2 * (
            (1.0 / 24.0) * (1.0 + 3.0 * cosSq2theta) * J(wr, tau, oneMinusS2) +
            (1.0 / 48.0) * (1.0 + 3.0 * cosSq2theta) * J(2.0 * wr, tau, oneMinusS2) +
            (3.0 / 4.0) * sinFourthTheta * (
                J(2.0 * we + wr, tau, oneMinusS2) +
                0.5 * J(2.0 * we + 2.0 * wr, tau, oneMinusS2) +
                0.5 * J(2.0 * we - 2.0 * wr, tau, oneMinusS2) +
                J(2.0 * we - wr, tau, oneMinusS2)
            ) +
            (3.0 / 8.0) * sinSq2theta * (
                J(we + wr, tau, oneMinusS2) +
                0.5 * J(we + 2.0 * wr, tau, oneMinusS2) +
                0.5 * J(we - 2.0 * wr, tau, oneMinusS2) +
                J(we - wr, tau, oneMinusS2)
            ) +
            (3.0 / 4.0) * (5.0 - cos2theta) * J(wI, tau, oneMinusS2) +
            1.5 * (3.0 + cos2theta) * J(2.0 * wI, tau, oneMinusS2)
        );
        return r1Rho;
    }
}
