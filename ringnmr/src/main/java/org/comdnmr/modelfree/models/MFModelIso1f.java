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
public class MFModelIso1f extends MFModelIso1 {

    double tauF;

    public MFModelIso1f(boolean fitTau, double targetTau, double tauFraction,
            boolean includeEx) {
        super(fitTau, targetTau, tauFraction, includeEx);
        nPars = includeEx ? 3 : 2;
    }

    public MFModelIso1f(double targetTau) {
        this(false, targetTau, 0.0, false);
    }

    public MFModelIso1f() {
        this(true, 0.0, 0.0, false);
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("Sf2", "Tau_f");
    }

    @Override
    public double[] calc(double[] omegas) {
        double[] J = new double[omegas.length];
        double sf2 = this.sf2 / sN;
        int j = 0;
        for (double omega : omegas) {
            omega *= 1.0e-9;
            double omega2 = omega * omega;
            double tauf = tauM * tauF / (tauM + tauF);
            double value1 = sf2 * tauM / (1.0 + omega2 * tauM * tauM);
            double value2 = (1.0 - sf2) * (tauf) / (1.0 + omega2 * tauf * tauf);
            J[j++] = 0.4e-9 * (value1 + value2);
        }
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
    }

    public double[] calc(double[] omegas, double s2, double tauF) {
        this.sf2 = s2;
        this.tauF = tauF;
        return calc(omegas);
    }

    @Override
    public double[] getStandardPars(double[] pars) {
        pars(pars);
        double[] stdPars = new double[5];
        stdPars[0] = tauM;
        stdPars[1] = sf2;
        stdPars[2] = tauF;
        stdPars[3] = 1.0;
        stdPars[4] = 0.0;
        return stdPars;
    }

    @Override
    public boolean checkParConstraints() {
        return tauF < tauM;
    }

    @Override
    public double[] getStart() {
        if (includeEx) {
            return getParValues(targetTau, 0.9, 0.015, 2.0);
        } else {
            return getParValues(targetTau, 0.9, 0.015);
        }
    }

    @Override
    public double[] getLower() {
        if (includeEx) {
            return getParValues(tauLower(), 0.0, 0.001, 0.0);
        } else {
            return getParValues(tauLower(), 0.0, 0.001);
        }
    }

    @Override
    public double[] getUpper() {
        if (includeEx) {
            return getParValues(tauUpper(), 1.0, SLOW_LIMIT, 100.0);
        } else {
            return getParValues(tauUpper(), 1.0, SLOW_LIMIT);

        }
    }

    @Override
    public int getNumber() {
        return 2;
    }

    @Override
    public String getName() {
        if (sN > 8.0) {
            return "modelD1f";
        } else {
            return "model1f";
        }
    }

}
