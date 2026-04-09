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
public class MFModelIso2s extends MFModelIso1s {

    double sf2;

    public MFModelIso2s(boolean fitTau, double targetTau, double tauFraction,
            boolean includeEx) {
        super(fitTau, targetTau, tauFraction, includeEx);
        nPars = includeEx ? 4 : 3;
    }

    public MFModelIso2s(double targetTau) {
        this(false, targetTau, 0.0, false);
    }

    public MFModelIso2s() {
        this(true, 0.0, 0.0, false);
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("Sf2", "Tau_s", "Ss2");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauM = 1.0e-9 * this.tauM;
        double tauS = 1.0e-9 * this.tauS;
        double ss2 = this.ss2;
        double sf2 = this.sf2 / sN;
        double s2 = ss2 * sf2;

        double tauMTimesPt4 = 0.4 * tauM;
        double tauM2 = tauM * tauM;
        double tauS2 = tauS * tauS;
        double tauM2TimesTauS2 = tauM2 * tauS2;
        double tauMPlusTauS = tauM + tauS;
        double tauMPlusTauS2 = tauMPlusTauS * tauMPlusTauS;

        double[] js = new double[omegas.length];
        int index = 0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double term1 = s2 / (1.0 + omega2 * tauM2);
            double term2 = sf2 * (1.0 - sf2) * (
                (tauS * tauMPlusTauS) / (omega2 * tauM2TimesTauS2 + tauMPlusTauS2)
            );
            js[index++] = tauMTimesPt4 * (term1 + term2);
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

        this.sf2 = pars[parStart];
        this.tauS = pars[parStart + 1];
        this.ss2 = pars[parStart + 2];
    }

    @Override
    public double[] getStandardPars(double[] pars) {
        pars(pars);
        return createStandardPars(sf2, 0.0, ss2, tauS);
    }

    public double[] calc(double[] omegas, double s2, double tauS, double sf2) {
        this.sf2 = s2;
        this.tauS = tauS;
        this.ss2 = sf2;
        return calc(omegas);
    }

    @Override
    public boolean checkParConstraints() {
        return tauS < tauM;
    }

    @Override
    public double[] getStart() {
        if (includeEx) {
            return getParValues(targetTau, 0.5, targetTau / 4.0, 0.5, 2.0);
        } else {
            return getParValues(targetTau, 0.5, targetTau / 4.0, 0.5);
        }
    }

    @Override
    public double[] getLower() {
        if (includeEx) {
            return getParValues(tauLower(), 0.0, SLOW_LIMIT, 0.0, 0.0);
        } else {
            return getParValues(tauLower(), 0.0, SLOW_LIMIT, 0.0);
        }
    }

    @Override
    public double[] getUpper() {
        if (includeEx) {
            return getParValues(tauUpper(), 1.0, targetTau / 2.0, 1.0, 100.0);
        } else {
            return getParValues(tauUpper(), 1.0, targetTau / 2.0, 1.0);

        }
    }

    @Override
    public int getNumber() {
        return 5;
    }

    @Override
    public String getName() {
        if (sN > 8.0) {
            return "modelD2s";
        } else {
            return "model2s";
        }
    }

}
