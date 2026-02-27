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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author brucejohnson
 */
public abstract class MFModelIso extends MFModel {
    double sN = 1.0;
    double tauM;
    double rEX;
    double tauFrac;

    public MFModelIso(boolean fitTau, double targetTau, double tauFraction,
                      boolean includeEx) {
        this.fitTau = fitTau;
        this.targetTau = targetTau;
        this.tauFrac = tauFraction;
        this.includeEx = includeEx;
        if (!fitTau) {
            tauM = targetTau;
        }
    }

    public MFModelIso(double targetTau) {
        this(false, targetTau, 0.0, false);
    }

    public MFModelIso() {
        this(true, 0.0, 0.0, false);
    }

    public List<String> getAllParNames(String... pars) {
        var parNames = new ArrayList<String>();
        if (fitTau) {
            parNames.add("Tau_e");
        }
        parNames.addAll(Arrays.asList(pars));
        if (includeEx) {
            parNames.add("Rex");
        }
        return parNames;
    }

    public void setSScale(double value) {
        sN = value;
    }

    public boolean fitTau() {
        return fitTau;
    }

    public double getTau() {
        return tauM;
    }

    public void setTauFraction(double value) {
        tauFrac = value;
    }

    public double tauLower() {
        return targetTau - targetTau * tauFrac;
    }

    public double tauUpper() {
        return targetTau + targetTau * tauFrac;
    }

    public abstract double[] getStart();

    public abstract double[] getStandardPars(double[] pars);


    public static MFModelIso buildModel(String modelName, boolean fitTau,
                                        double tau, double tauFrac,
                                        boolean fitExchange) {
        MFModelIso model;
        if (modelName.startsWith("model")) {
            modelName = modelName.substring(5);
        }
        model = switch (modelName) {
            case "1", "D1" -> new MFModelIso1(fitTau, tau, tauFrac, fitExchange);
            case "1f", "D1f" -> new MFModelIso1f(fitTau, tau, tauFrac, fitExchange);
            case "1s", "D1s" -> new MFModelIso1s(fitTau, tau, tauFrac, fitExchange);
            case "2s", "D2s" -> new MFModelIso2s(fitTau, tau, tauFrac, fitExchange);
            case "2f", "D2f" -> new MFModelIso2f(fitTau, tau, tauFrac, fitExchange);
            case "1sf", "D1sf" -> new MFModelIso1sf(fitTau, tau, tauFrac, fitExchange);
            case "2sf", "D2sf" -> new MFModelIso2sf(fitTau, tau, tauFrac, fitExchange);
            case "2sfx", "D2sfx" -> new MFModelIso2sf(fitTau, tau, tauFrac, fitExchange);
            default -> throw new IllegalArgumentException("Unknown model " + modelName);
        };
        if (modelName.charAt(0) == 'D') {
            model.setSScale(9.0);
        }
        return model;

    }

}
