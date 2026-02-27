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

import org.comdnmr.modelfree.RelaxEquations;
import org.comdnmr.util.CoMDPreferences;

/**
 *
 * @author Simon Hulse
 */
public enum SSR1RhoEquation implements EquationType {

    CSA("CSA", 0, "tauc", "s2") {

        @Override
        // nu1 and nuR are both provided in kHz
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            // FIXME: During simulations, getting map = [1, 0]... not sure why
            // expected tauc = par[map[0]] and s2 = par[map[1]]
            double tauc = par[map[0]];
            double s2 = par[map[1]];
            double nu1kHz = x[0];
            double nuRkHz = x[3];
            double omega1 = 2.0e3 * Math.PI * nu1kHz;
            double omegaR = 2.0e3 * Math.PI * nuRkHz;

            double b0 = 1.0e6 * CoMDPreferences.getRefField();
            // TODO: currently fixed to carbon
            RelaxEquations relaxEquations = new RelaxEquations(b0, "H", "C");
            double r1rho = relaxEquations.r1RhoCSA(omegaR, omega1, tauc, s2);
            return r1rho;
        }
    },

    DIPOLAR_IS("DIPOLAR_IS", 0, "tauc", "s2") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            double tauc = par[map[0]];
            double s2 = par[map[1]];
            double nu1kHz = x[0];
            double nuRkHz = x[3];
            double omega1 = 2.0e3 * Math.PI * nu1kHz;
            double omegaR = 2.0e3 * Math.PI * nuRkHz;

            double b0 = 1.0e6 * CoMDPreferences.getRefField();
            // TODO: Get nuclei
            RelaxEquations relaxEquations = new RelaxEquations(b0, "H", "C");
            double r1rho = relaxEquations.r1RhoIS(omegaR, omega1, tauc, s2);
            return r1rho;
        }
    },

    DIPOLAR_AB("DIPOLAR_AB", 0, "tauc", "s2") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            double tauc = par[map[0]];
            double s2 = par[map[1]];
            double nu1kHz = x[0];
            double nuRkHz = x[3];
            double omega1 = 2.0e3 * Math.PI * nu1kHz;
            double omegaR = 2.0e3 * Math.PI * nuRkHz;

            double b0 = 1.0e6 * CoMDPreferences.getRefField();
            // TODO: Get nuclei
            RelaxEquations relaxEquations = new RelaxEquations(b0, "H", "H");
            double r1rho = relaxEquations.r1RhoAB(omegaR, omega1, tauc, s2);
            return r1rho;
        }
    },

    DIPOLAR_AA("DIPOLAR_AA", 0, "tauc", "s2") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            double tauc = par[map[0]];
            double s2 = par[map[1]];
            double nu1kHz = x[0];
            double nuRkHz = x[3];
            double omega1 = 2.0e3 * Math.PI * nu1kHz;
            double omegaR = 2.0e3 * Math.PI * nuRkHz;

            double b0 = 1.0e6 * CoMDPreferences.getRefField();
            // TODO: Get nuclei
            RelaxEquations relaxEquations = new RelaxEquations(b0, "H", "H");
            double r1rho = relaxEquations.r1RhoAA(omegaR, omega1, tauc, s2);
            return r1rho;
        }
    };

    final String equationName;
    final int nGroupPars;
    String[] parNames;
    static final String[] equationNames = {"CSA", "DIPOLAR_IS", "DIPOLAR_AB", "DIPOLAR_AA"};

    SSR1RhoEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    public static String[] getAllEquationNames() {
        return equationNames;
    }

    @Override
    public String getName() {
        return equationName;
    }

    @Override
    public String[] getParNames() {
        return parNames;
    }

    @Override
    public int getNGroupPars() {
        return nGroupPars;
    }

    public double getMinX() {
        return 1.0;
    }

    public double getMaxX() {
        return 50.0;
    }

    @Override
    public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
        // TODO: should this be X.length or X[0].length?
        int n = X[0].length;
        double[] yCalc = new double[n];
        double[] x = new double[X.length];
        for (int i = 0; i < n; i++) {
            for (int j=0;j<x.length;j++) {
                x[j] = X[j][i];
            }
            yCalc[i] = calculate(par, map, x, idNum);
        }
        return yCalc;
    }

    // TODO: integrate with `map`
    public double[][] boundaries(
        double[] guesses,
        double[][] xValues,
        double[] yValues,
        int[][] map,
        int[] idNums,
        int nID) {
        return new double[][]{{0.0, 0.0}, {1.0e-3, 1.0}};
    }

    @Override
    public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
        return new double[]{1.0e-6, 1.0};
    }

    @Override
    public double getRex(double[] pars, int[] map, double field) {
        return 0.0;
    }

    @Override
    public double getKex(double[] pars) {
        return pars[0];
    }

    @Override
    public double getKex(double[] pars, int id) {
        return pars[0];
    }

    // TODO: For Bruce
    @Override
    public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
        return makeMap(1);
    }
}
