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
    double complexityS = 0.0;
    double complexityTau = 0.0;

    public MFModelIso1sf(boolean fitTau, double targetTau, double tauFraction,
                         boolean includeEx) {
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
        double[] J = new double[omegas.length];
        int j = 0;
        double ss2 = this.ss2;
        double sf2 = 1.0 / sN;
        double s2 = ss2 * sf2;
        for (double omega : omegas) {
            omega *= 1.0e-9;
            double omega2 = omega * omega;
            double vM = s2 * tauM / (1.0 + omega2 * tauM * tauM);
            double vMS = sf2 * (1.0 - ss2) * (tauM * tauS * (tauM + tauS))
                    / (tauM * tauM * tauS * tauS * omega2 + (tauM + tauS) * (tauM + tauS));
            double vMF = (1.0 - sf2) * ss2 * (tauM * tauF * (tauM + tauF))
                    / (tauM * tauM * tauF * tauF * omega2 + (tauM + tauF) * (tauM + tauF));
            double tauMFS = tauF * (tauM + tauS) + tauM * tauS;
            double vMFS = (1.0 - sf2) * (1.0 - ss2) * (tauF * tauM * tauS * tauMFS)
                    / (tauF * tauF * tauM * tauM * tauS * tauS * omega2
                    + tauMFS * tauMFS);
            J[j++] = 0.4e-9 * (vM + vMS + vMF + vMFS);
        }
        complexityS
                = Math.abs(1.0 - ss2);
        complexityTau
                = (Math.log10(tauS + 0.001) + 3.0)
                + (Math.log10(tauF + 0.001) + 3.0);
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
            return getParValues(targetTau, 0.1, 0.9, 0.01, 2.0);
        } else {
            return getParValues(targetTau, SLOW_LIMIT / 5.0, 0.9, SLOW_LIMIT * 5.0);
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
