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
 * @author brucejohnson
 */
public class MFModelIso1sf extends MFModelIso2f {
    double tauS;

    public MFModelIso1sf(
        boolean fitTau,
        double targetTau,
        double tauFraction,
        boolean includeEx
    ) {
        super(fitTau, targetTau, tauFraction, includeEx);
        nPars = includeEx ? 4 : 3;
    }

    public MFModelIso1sf(double targetTau) {
        this(false, targetTau, 0.0, false);
    }

    public MFModelIso1sf() {
        this(true, 0.0, 0.0, false);
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("Tau_f", "Ss2", "Tau_s");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauM = 1.0e-9 * this.tauM;
        double tauF = 1.0e-9 * this.tauF;
        double tauS = 1.0e-9 * this.tauS;
        double ss2 = this.ss2;
        double sf2 = 1.0 / sN;
        double s2 = ss2 * sf2;

        double tauMTimesPt4 = 0.4 * tauM;
        double tauM2 = tauM * tauM;
        double tauF2 = tauF * tauF;
        double tauS2 = tauS * tauS;
        double tauMPlusTauF = tauM + tauF;
        double tauMPlusTauF2 = tauMPlusTauF * tauMPlusTauF;
        double tauMPlusTauS = tauM + tauS;
        double tauMPlusTauS2 = tauMPlusTauS * tauMPlusTauS;
        double tauFTimesTauS = tauF * tauS;
        double tauPrime = tauF * (tauM + tauS) + tauM * tauS;
        double tauPrime2 = tauPrime * tauPrime;
        double tauM2TimesTauF2 = tauM2 * tauF2;
        double tauM2TimesTauS2 = tauM2 * tauS2;
        double tauM2TimesTauF2TimesTauS2 = tauM2 * tauF2 * tauS2;

        double[] js = new double[omegas.length];
        int index = 0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double term1 = s2 / (1.0 + omega2 * tauM2);
            double term2 = sf2 * (1.0 - ss2) * (
                (tauS * tauMPlusTauS) / (omega2 * tauM2TimesTauS2 + tauMPlusTauS2)
            );
            double term3 = (1.0 - sf2) * ss2 * (
                (tauF * tauMPlusTauF) / (omega2 * tauM2TimesTauF2 + tauMPlusTauF2)
            );
            double term4 = (1.0 - sf2) * (1.0 - ss2) * (
                (tauFTimesTauS * tauPrime) / (omega2 * tauM2TimesTauF2TimesTauS2 + tauPrime2)
            );
            js[index++] = tauMTimesPt4 * (term1 + term2 + term3 + term4);
        }
        return js;
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

        this.tauF = pars[parStart];
        this.ss2 = pars[parStart + 1];
        this.tauS = pars[parStart + 2];
    }

    public double[] sortPars(double[] pars) {
        pars(pars);
        double[] sortPars = new double[4];
        sortPars[0] = tauM;
        boolean swapIt = tauF > tauS;

        if (swapIt) {
            sortPars[1] = tauS;
            sortPars[2] = ss2;
            sortPars[3] = tauF;
        } else {
            sortPars[1] = tauF;
            sortPars[2] = ss2;
            sortPars[3] = tauS;
        }
        return sortPars;
    }

    @Override
    public double[] getStandardPars(double[] pars) {
        double[] sortPars = sortPars(pars);
        double tauF = sortPars[1];
        double ss2 = sortPars[2];
        double tauS = sortPars[3];
        return createStandardPars(1.0, tauF, ss2, tauS);
    }

    @Override
    public double[] calc(double[] omegas, double tauF, double ss2, double tauS) {
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
            return getParValues(targetTau, SLOW_LIMIT / 2.0, 0.5, targetTau / 4.0, 2.0);
        } else {
            return getParValues(targetTau, SLOW_LIMIT / 2.0, 0.5, targetTau / 4.0);
        }
    }

    @Override
    public double[] getLower() {
        if (includeEx) {
            return getParValues(tauLower(), 0.001, 0.0, SLOW_LIMIT, 0.0);
        } else {
            return getParValues(tauLower(), 0.001, 0.0, SLOW_LIMIT);
        }
    }

    @Override
    public double[] getUpper() {
        if (includeEx) {
            return getParValues(tauUpper(), SLOW_LIMIT, 1.0, targetTau / 2.0, 100.0);
        } else {
            return getParValues(tauUpper(), SLOW_LIMIT, 1.0, targetTau / 2.0);

        }
    }

    @Override
    public int getNumber() {
        return 6;
    }

    @Override
    public String getName() {
        if (sN > 8.0) {
            return "modelD1sf";
        } else {
            return "model1sf";
        }
    }

}
