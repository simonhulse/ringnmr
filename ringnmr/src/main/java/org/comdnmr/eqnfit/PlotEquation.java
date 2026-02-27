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

import org.comdnmr.fit.ResidueFitter;
import org.comdnmr.modelfree.RelaxEquations;
import org.comdnmr.modelfree.models.MFModelIso;

/**
 * @author Bruce Johnson
 */
public class PlotEquation {

    String expType;
    String name;
    double[] pars;
    double[] errs;
    double[] extras;

    public PlotEquation(String expType, String name, double[] pars, double[] errs, double[] extras) {
        this.expType = expType;
        this.name = name;
        this.pars = pars.clone();
        this.errs = errs.clone();
        this.extras = extras.clone();
    }

    @Override
    public PlotEquation clone() {
        return new PlotEquation(expType, name, pars, errs, extras);
    }

    public String getName() {
        return name;
    }

    public void setExtra(double[] extras) {
        this.extras = extras.clone();
    }

    public double getExtra(int i) {
        return extras[i];
    }

    public double[] getExtras() {
        return extras;
    }

    public double[] getPars() {
        return pars;
    }

    public double[] getErrs() {
        return errs;
    }

    public double calculate(double[] simPars, double[] xValue) {
        double val;
        if (expType.startsWith("model")) {
            val = calculateSpectralDensity(xValue, simPars);
        } else {
            EquationType equationType = ResidueFitter.getEquationType(expType, name);
            int[][] map = equationType.makeMap(1);
            val = equationType.calculate(simPars, map[0], xValue, 0);
            if (expType.equals("ssr1rho")) {
                val = Math.log10(val);
            }
        }
        return val;
    }

    public double calculate(double[] xValue) { return calculate(pars, xValue); }

    private double calculateSpectralDensity(double[] xValue) {
        return calculateSpectralDensity(xValue, pars);
    }

    private double calculateSpectralDensity(double[] xValue, double[] simPars) {
        var model = MFModelIso.buildModel(expType, true, 0.0, 0.0, false);
        double[] omegas = {xValue[0] * 1.0e9};
        double[] specDens = model.calc(omegas, simPars);
        return Math.log10(specDens[0] * 1.0e9);
    }

    private double calculateR(double[] xValue) {
        var model = MFModelIso.buildModel(expType, true, 0.0, 0.0, false);
        RelaxEquations relaxEquations = new RelaxEquations(800.0e6, "H", "N");

        double[] omegas = relaxEquations.getOmegas();
        double[] specDens = model.calc(omegas, pars);
        return relaxEquations.R1(specDens);
    }

    public double calculate(double xValue) {
        EquationType equationType = ResidueFitter.getEquationType(expType, name);
        int[][] map = equationType.makeMap(1);
        double[] ax = new double[extras.length];
        ax[0] = xValue;
        System.arraycopy(extras, 1, ax, 1, extras.length - 1);
        return calculate(ax);
    }

    public double getMinX() {
        if (expType.startsWith("model")) {
            return 0.0;
        } else {
            EquationType equationType = ResidueFitter.getEquationType(expType, name);
            if (equationType == null) {
                return 0.0;
            }
            return equationType.getMinX();
        }
    }

    public double getMaxX() {
        if (expType.startsWith("model")) {
            if (expType.contains("D")) {
                return 2.0;
            } else {
                return 5.0;
            }
        } else {
            EquationType equationType = ResidueFitter.getEquationType(expType, name);
            if (equationType == null) {
                return 0.0;
            }
            return equationType.getMaxX();
        }
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append(name);
        for (double par : pars) {
            builder.append(" ").append(par);
        }
        for (double extra : extras) {
            builder.append(" ").append(extra);
        }
        return builder.toString();

    }

}
