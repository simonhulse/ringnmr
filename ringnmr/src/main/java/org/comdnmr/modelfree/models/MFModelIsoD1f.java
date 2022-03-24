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
public class MFModelIsoD1f extends MFModelIso1 {

    double tauF;

    public MFModelIsoD1f(boolean fitTau, double targetTau, double tauFraction,
                         boolean includeEx) {
        super(fitTau, targetTau, tauFraction, includeEx);
        nPars = includeEx ? 3 : 2;
    }

    public MFModelIsoD1f(double targetTau) {
        this(false, targetTau, 0.0, false);
    }

    public MFModelIsoD1f() {
        this(true, 0.0, 0.0, false);
    }

    @Override
    public List<String> getParNames() {
        return getAllParNames("Sf2", "Tau_f");
    }

    @Override
    public double[] calc(double[] omegas) {
        double tauMx = tauM * 1.0e-9;
        double tauFx = tauF * 1.0e-9;
        double sf2p = sf2 / 9.0;
        double[] J = new double[omegas.length];
        int j = 0;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double tauf = tauMx * tauFx / (tauMx + tauFx);
            double value1 = sf2p * tauMx / (1.0 + omega2 * tauMx * tauMx);
            double value2 = (1.0 - sf2p) * (tauf) / (1.0 + omega2 * tauf * tauf);
            J[j++] = 0.4 * (value1 + value2);
        }
        return J;
    }

    @Override
    public double[] calc(double[] omegas, double[] pars) {
        int parStart = 0;
        if (fitTau) {
            tauM = pars[0];
            parStart = 1;
        }

        this.sf2 = pars[parStart];
        this.tauF = pars[parStart + 1];
        return calc(omegas);
    }

    public double[] calc(double[] omegas, double s2, double tauF) {
        this.sf2 = s2;
        this.tauF = tauF;
        return calc(omegas);
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
        return "model1f";
    }

}