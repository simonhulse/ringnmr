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
package org.comdnmr.modelfree.models;

import java.util.List;

/**
 *
 * @author brucejohnson
 */
public class MFModelIso2sf extends MFModelIso2f {
    protected static double tauPrime = 30.0e-3;
    double tauS;
    double complexityS = 0.0;
    double complexityTau = 0.0;

    public MFModelIso2sf(boolean fitTau, double targetTau, double tauFraction,
            boolean includeEx) {
        super(fitTau, targetTau, tauFraction, includeEx);
        nPars = includeEx ? 5 : 4;
    }

    public MFModelIso2sf(double targetTau) {
        this(false, targetTau, 0.0, false);
    }

    public MFModelIso2sf() {
        this(true, 0.0, 0.0, false);
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("Sf2", "Tau_f", "Ss2", "Tau_s");
    }

    @Override
    public double[] calc(double[] omegas) {
        double[] J = new double[omegas.length];
        int j = 0;
        double ss2 = this.ss2;
        double sf2 = this.sf2 / sN;
        double s2 = ss2 * sf2;
        for (double omega : omegas) {
            omega *= 1.0e-9;
            double omega2 = omega * omega;
            double vM = s2 / (1.0 + omega2 * tauM * tauM);
            double vMF = ((1.0 - sf2) * tauF * (tauM + tauF))
                    / (omega2 * tauM * tauM * tauF * tauF + (tauM + tauF) * (tauM + tauF));
            double vMS = (sf2 * (1.0 - ss2) * tauS * (tauM + tauS))
                    / (omega2 * tauM * tauM * tauS * tauS + (tauM + tauS) * (tauM + tauS));
            J[j++] = 0.4 * tauM * 1.0e-9 *  (vM + vMF + vMS);
        }

        complexityS = Math.abs(1.0 - sf2) + Math.abs(1.0 - ss2);
        complexityTau =
                Math.log10((tauS + tauPrime) / tauPrime) +
                        Math.log10((tauF + tauPrime) / tauPrime);
        return J;
    }

    @Override
    public double[] calc(double[] omegas, double[] pars) {
        pars(pars);
        return calc(omegas);
    }

    @Override
    public void pars(double[] pars) {
        int parStart = 0;
        if (fitTau) {
            tauM = pars[0];
            parStart = 1;
        }

        this.sf2 = pars[parStart];
        this.tauF = pars[parStart + 1];
        this.ss2 = pars[parStart + 2];
        this.tauS = pars[parStart + 3];
    }

    public double[] sortPars(double[] pars) {
        pars(pars);
        double[] sortPars = new double[5];
        sortPars[0] = tauM;
        double tauLimit = 7.0e-3;
        boolean swapIt = false;
        if ((tauS < tauLimit) && (tauF < tauLimit)) {
            if (ss2 < sf2) {
                swapIt = true;
            }
        } else if (tauF > tauS) {
            swapIt = true;
        }

        if (swapIt) {
            sortPars[1] = ss2;
            sortPars[2] = tauS;
            sortPars[3] = sf2;
            sortPars[4] = tauF;
        } else {
            sortPars[1] = sf2;
            sortPars[2] = tauF;
            sortPars[3] = ss2;
            sortPars[4] = tauS;
        }
        return sortPars;
    }

    @Override
    public double[] getStandardPars(double[] pars) {
        return sortPars(pars);
    }

    @Override
    public double getComplexityS() {
        return complexityS;
    }

    @Override
    public double getComplexityTau() {
        return complexityTau;
    }

    public double[] calc(double[] omegas, double sf2, double tauF, double ss2, double tauS) {
        this.sf2 = sf2;
        this.tauF = tauF;
        this.ss2 = ss2;
        this.tauS = tauS;
        return calc(omegas);
    }

    @Override
    public boolean checkParConstraints() {
        return tauF < tauM && tauS < tauM;
    }

    @Override
    public double[] getStart() {
        if (includeEx) {
            return getParValues(targetTau, 0.9, 0.1, 0.9, 0.01, 2.0);
        } else {
            return getParValues(targetTau, 0.9, SLOW_LIMIT / 5.0, 0.9, SLOW_LIMIT * 5.0);
        }
    }

    @Override
    public double[] getLower() {
        if (includeEx) {
            return getParValues(tauLower(), 0.0, 0.000001, 0.0, SLOW_LIMIT, 0.0);
        } else {
            return getParValues(tauLower(), 0.0, 0.000001, 0.0, SLOW_LIMIT);
        }
    }

    @Override
    public double[] getUpper() {
        if (includeEx) {
            return getParValues(tauUpper(), 1.0, SLOW_LIMIT, 1.0, targetTau / 2.0, 100.0);
        } else {
            return getParValues(tauUpper(), 1.0, SLOW_LIMIT, 1.0, targetTau / 2.0);

        }
    }

    @Override
    public int getNumber() {
        return 7;
    }

    @Override
    public String getName() {
        if (sN > 8.0) {
            return "modelD2sf";
        } else {
            return "model2sf";
        }
    }

}
