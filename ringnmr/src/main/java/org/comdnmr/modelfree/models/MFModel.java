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
public abstract class MFModel {

    public final static double SLOW_LIMIT = 0.15;
    int nPars;
    boolean fitTau;
    boolean includeEx;
    double targetTau;

    public abstract double[] calc(double[] omegas, double[] pars);

    public abstract double[] calc(double[] omega);

    public double[] getParValues(double... parValues) {
        int n = getNPars();
        int start = fitTau ? 0 : 1;
        double[] values = new double[n];
        System.arraycopy(parValues, start, values, 0, n);
        return values;
    }

    public abstract List<String> getParNames();

    public abstract double[] getLower();

    public abstract double[] getUpper();

    public boolean checkParConstraints() {
        return true;
    }

    public double calcWeight() {
        return 0.0;
    }

    public double getComplexityS2F() {
        return 0.0;
    }

    public double getComplexityS2S() {
        return 0.0;
    }

    public double getComplexityTauF() {
        return 0.0;
    }

    public double getComplexityTauS() {
        return 0.0;
    }

    public boolean includesEx() {
        return includeEx;
    }

    public int getNPars() {
        return (fitTau) ? nPars + 1 : nPars;
    }

    public abstract int getNumber();

    public String getName() {
        return "model" + getNumber();
    }
}
