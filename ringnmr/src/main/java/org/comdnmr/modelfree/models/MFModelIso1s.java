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
public class MFModelIso1s extends MFModelIso {

    double ss2;
    double tauS;

    public MFModelIso1s(boolean fitTau, double targetTau, double tauFraction,
            boolean includeEx) {
        super(fitTau, targetTau, tauFraction, includeEx);
        nPars = includeEx ? 3 : 2;
    }

    public MFModelIso1s(double targetTau) {
        this(false, targetTau, 0.0, false);
    }

    public MFModelIso1s() {
        this(true, 0.0, 0.0, false);
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("Ss2", "Tau_s");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauM = 1.0e-9 * this.tauM;
        double tauS = 1.0e-9 * this.tauS;
        double ss2 = this.ss2;

        double tauMTimesPt4 = 0.4 * tauM;
        double tauM2 = tauM * tauM;
        double tauS2 = tauS * tauS;
        double tauM2TimesTauS2 = tauM2 * tauS2;
        double tauMPlusTauS = tauM + tauS;
        double tauMPlusTauS2 = tauMPlusTauS * tauMPlusTauS;

        int index = 0;
        double[] js = new double[omegas.length];
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double term1 = ss2 / (1.0 + omega2 * tauM2);
            double term2 = ((1.0 - ss2) * tauS * tauMPlusTauS) /
                (omega2 * tauM2TimesTauS2 + tauMPlusTauS2);
            js[index++] = tauMTimesPt4 * (term1 + term2);
        }
        return js;
    }

    @Override
    public double[] calc(double[] omegas, double[] pars) {
        pars(pars);
        return calc(omegas);
    }

    public void pars(double[] pars) {
        int parStart = 0;
        if (fitTau) {
            tauM = pars[0];
            parStart = 1;
        }

        this.ss2 = pars[parStart];
        this.tauS = pars[parStart + 1];
    }

    @Override
    public double[] getStandardPars(double[] pars) {
        pars(pars);
        double[] stdPars = new double[5];
        stdPars[0] = tauM;
        stdPars[1] = 1.0;
        stdPars[2] = 0.0;
        stdPars[3] = ss2;
        stdPars[4] = tauS;
        return stdPars;
    }

    public double[] calc(double[] omegas, double s2, double tauS) {
        this.ss2 = s2;
        this.tauS = tauS;
        return calc(omegas);
    }

    @Override
    public boolean checkParConstraints() {
        return tauS < tauM;
    }

    @Override
    public double[] getStart() {
        if (includeEx) {
            return getParValues(targetTau, 0.5, targetTau / 4.0, 2.0);
        } else {
            return getParValues(targetTau, 0.5, targetTau / 4.0);
        }
    }

    @Override
    public double[] getLower() {
        if (includeEx) {
            return getParValues(tauLower(), 0.0, SLOW_LIMIT, 0.0);
        } else {
            return getParValues(tauLower(), 0.0, SLOW_LIMIT);
        }
    }

    @Override
    public double[] getUpper() {
        if (includeEx) {
            return getParValues(tauUpper(), 1.0, targetTau / 2.0, 100.0);
        } else {
            return getParValues(tauUpper(), 1.0, targetTau / 2.0);

        }
    }

    @Override
    public int getNumber() {
        return 3;
    }

    @Override
    public String getName() {
        if (sN > 8.0) {
            return "modelD1s";
        } else {
            return "model1s";
        }
    }
}
