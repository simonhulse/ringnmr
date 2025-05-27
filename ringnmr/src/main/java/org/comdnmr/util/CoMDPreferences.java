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
package org.comdnmr.util;

import org.comdnmr.data.Experiment;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.prefs.Preferences;

/**
 *
 * @author brucejohnson
 */
public class CoMDPreferences {

    static private Double refField = null;
    static private Double cpmgMaxFreq = null;
    static private Double rexRatio = null;
    static private Double startRadius = null;
    static private Double finalRadius = null;
    static private Double tolerance = null;
    static private Boolean weightFit = null;
    static private Boolean absValueFit = null;
    static private Boolean nonParametricBootstrap = null;
    static private Boolean neuralNetworkGuess = null;
    static private Boolean calR1rhoCorr = null;
    static private Integer sampleSize = null;
    static private String optimizer = null;
    static private String bootStrapOptimizer = null;
    private static Map<String, Boolean> cestEqnMap = null;
    private static Map<String, Boolean> cpmgEqnMap = null;
    private static Map<String, Boolean> r1rhoEqnMap = null;
    private static Map<String, Boolean> expEqnMap = null;
    private static Map<String, Boolean> noeEqnMap = null;
    private static final String DEFAULT_CPMG_EQNS = "NOEX;true\nCPMGFAST;true\nCPMGSLOW;true";
    private static final String DEFAULT_CEST_EQNS = "NOEX;true\nTROTT_PALMER;true";
    private static final String DEFAULT_R1RHO_EQNS = "NOEX;true\nTROTT_PALMER;true";
    private static final String DEFAULT_EXP_EQNS = "EXPAB;true";
    private static final String DEFAULT_NOE_EQNS = "NOE;true";
    private static Double deltaABdiff = null;

    static Preferences getPrefs() {
        Preferences prefs = Preferences.userNodeForPackage(CoMDPreferences.class);
        return prefs;
    }

    public static Double getCPMGMaxFreq() {
        if (cpmgMaxFreq == null) {
            String value = getPrefs().get("CPMG_MAX_FREQ", "2000.0");
            cpmgMaxFreq = Double.parseDouble(value);
        }
        return cpmgMaxFreq;
    }

    public static void setCPMGMaxFreq(Double value) {
        cpmgMaxFreq = value;
        if (value != null) {
            getPrefs().put("CPMG_MAX_FREQ", value.toString());
        } else {
            getPrefs().remove("CPMG_MAX_FREQ");
        }
    }

    public static void setRefField(Double value) {
        refField = value;
        if (value != null) {
            getPrefs().put("REF_FIELD", value.toString());
        } else {
            getPrefs().remove("REF_FIELD");
        }
    }

    public static Double getRefField() {
        if (refField == null) {
            String value = getPrefs().get("REF_FIELD", String.valueOf(CoMDDefaults.getRefField()));
            refField = Double.parseDouble(value);
        }
        return refField;
    }

    public static Integer getSampleSize() {
        if (sampleSize == null) {
            String value = getPrefs().get("SAMPLE_SIZE", String.valueOf(CoMDDefaults.getSampleSize()));
            sampleSize = Integer.parseInt(value);
        }
        return sampleSize;
    }

    public static void setSampleSize(Integer value) {
        sampleSize = value;
        if (value != null) {
            getPrefs().put("SAMPLE_SIZE", value.toString());
        } else {
            getPrefs().remove("SAMPLE_SIZE");
        }
    }

    public static String getOptimizer() {
        if (optimizer == null) {
            String value = getPrefs().get("OPTIMIZER", CoMDDefaults.getOptimizer());
            optimizer = value;
        }
        return optimizer;
    }

    public static void setOptimizer(String value) {
        optimizer = value;
        if (value != null) {
            getPrefs().put("OPTIMIZER", value);
        } else {
            getPrefs().remove("OPTIMIZER");
        }
    }

    public static String getBootStrapOptimizer() {
        if (bootStrapOptimizer == null) {
            String value = getPrefs().get("BOOTSTRAP-OPTIMIZER", CoMDDefaults.getOptimizer());
            bootStrapOptimizer = value;
        }
        return bootStrapOptimizer;
    }

    public static void setBootStrapOptimizer(String value) {
        bootStrapOptimizer = value;
        if (value != null) {
            getPrefs().put("BOOTSTRAP-OPTIMIZER", value);
        } else {
            getPrefs().remove("BOOTSTRAP-OPTIMIZER");
        }
    }

    public static Double getRexRatio() {
        if (rexRatio == null) {
            String value = getPrefs().get("REX_RATIO", String.valueOf(CoMDDefaults.getRexRatio()));
            rexRatio = Double.parseDouble(value);
        }
        return rexRatio;
    }

    public static void setRexRatio(Double value) {
        rexRatio = value;
        if (value != null) {
            getPrefs().put("REX_RATIO", value.toString());
        } else {
            getPrefs().remove("REX_RATIO");
        }
    }

    public static Double getDeltaABDiff() {
        if (deltaABdiff == null) {
            String value = getPrefs().get("DELTAAB_DIFF", String.valueOf(CoMDDefaults.getDeltaABDiff()));
            deltaABdiff = Double.parseDouble(value);
        }
        return deltaABdiff;
    }

    public static void setDeltaABDiff(Double value) {
        deltaABdiff = value;
        if (value != null) {
            getPrefs().put("DELTAAB_DIFF", value.toString());
        } else {
            getPrefs().remove("DELTAAB_DIFF");
        }
    }

    public static Double getStartRadius() {
        if (startRadius == null) {
            String value = getPrefs().get("START_RADIUS", String.valueOf(CoMDDefaults.getStartRadius()));
            startRadius = Double.parseDouble(value);
        }
        return startRadius;
    }

    public static void setStartingRadius(Double value) {
        startRadius = value;
        if (value != null) {
            getPrefs().put("START_RADIUS", value.toString());
        } else {
            getPrefs().remove("START_RADIUS");
        }
    }

    public static Double getFinalRadius() {
        if (finalRadius == null) {
            String value = getPrefs().get("FINAL_RADIUS", String.valueOf(CoMDDefaults.getFinalRadius()));
            finalRadius = Double.parseDouble(value);
        }
        return finalRadius;
    }

    public static void setFinalRadius(Double value) {
        finalRadius = value;
        if (value != null) {
            getPrefs().put("FINAL_RADIUS", value.toString());
        } else {
            getPrefs().remove("FINAL_RADIUS");
        }
    }

    public static Double getTolerance() {
        if (tolerance == null) {
            String value = getPrefs().get("TOLERANCE", String.valueOf(CoMDDefaults.getTolerance()));
            tolerance = Double.parseDouble(value);
        }
        return tolerance;
    }

    public static void setTolerance(Double value) {
        tolerance = value;
        if (value != null) {
            getPrefs().put("TOLERANCE", value.toString());
        } else {
            getPrefs().remove("TOLERANCE");
        }
    }

    public static Boolean getWeightFit() {
        if (weightFit == null) {
            String value = getPrefs().get("WEIGHT_FIT", String.valueOf(CoMDDefaults.getWeightFit()));
            weightFit = Boolean.parseBoolean(value);
        }
        return weightFit;
    }

    public static void setWeightFit(Boolean value) {
        weightFit = value;
        if (value != null) {
            getPrefs().put("WEIGHT_FIT", value.toString());
        } else {
            getPrefs().remove("WEIGHT_FIT");
        }
    }

    public static Boolean getAbsValueFit() {
        if (absValueFit == null) {
            String value = getPrefs().get("ABSVAL_FIT", String.valueOf(CoMDDefaults.getAbsValueFit()));
            absValueFit = Boolean.parseBoolean(value);
        }
        return absValueFit;
    }

    public static void setAbsValueFit(Boolean value) {
        absValueFit = value;
        if (value != null) {
            getPrefs().put("ABSVAL_FIT", value.toString());
        } else {
            getPrefs().remove("ABSVAL_FIT");
        }
    }

    public static Boolean getNonParametricBootstrap() {
        if (nonParametricBootstrap == null) {
            String value = getPrefs().get("NONPARAMETRIC", String.valueOf(CoMDDefaults.getNonParametricBootstrap()));
            nonParametricBootstrap = Boolean.parseBoolean(value);
        }
        return nonParametricBootstrap;
    }

    public static void setNonParametric(Boolean value) {
        nonParametricBootstrap = value;
        if (value != null) {
            getPrefs().put("NONPARAMETRIC", value.toString());
        } else {
            getPrefs().remove("NONPARAMETRIC");
        }
    }

    public static Boolean getNeuralNetworkGuess() {
        if (neuralNetworkGuess == null) {
            String value = getPrefs().get("NEURALNETWORKGUESS", String.valueOf(CoMDDefaults.getNeuralNetworkGuess()));
            neuralNetworkGuess = Boolean.parseBoolean(value);
        }
        return neuralNetworkGuess;
    }

    public static void setNeuralNetworkGuess(Boolean value) {
        neuralNetworkGuess = value;
        if (value != null) {
            getPrefs().put("NEURALNETWORKGUESS", value.toString());
        } else {
            getPrefs().remove("NEURALNETWORKGUESS");
        }
    }

    public static Boolean getCalR1rhoCorr() {
        if (calR1rhoCorr == null) {
            String value = getPrefs().get("CALR1RHOCORR", String.valueOf(CoMDDefaults.getCalR1rhoCorr()));
            calR1rhoCorr = Boolean.parseBoolean(value);
        }
        return calR1rhoCorr;
    }

    public static void setCalR1rhoCorr(Boolean value) {
        calR1rhoCorr = value;
        if (value != null) {
            getPrefs().put("CALR1RHOCORR", value.toString());
        } else {
            getPrefs().remove("CALR1RHOCORR");
        }
    }

    static Map<String, Boolean> getCPMGEqnMap() {
        if (cpmgEqnMap == null) {
            Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
            String eqns = prefs.get("CPMG_EQNS", DEFAULT_CPMG_EQNS);
            cpmgEqnMap = stringToMap(eqns);
        }
        return cpmgEqnMap;
    }

    static Map<String, Boolean> getCESTEqnMap() {
        if (cestEqnMap == null) {
            Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
            String eqns = prefs.get("CEST_EQNS", DEFAULT_CEST_EQNS);
            cestEqnMap = stringToMap(eqns);
        }
        return cestEqnMap;
    }

    static Map<String, Boolean> getR1RhoEqnMap() {
        if (r1rhoEqnMap == null) {
            Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
            String eqns = prefs.get("R1RHO_EQNS", DEFAULT_R1RHO_EQNS);
            r1rhoEqnMap = stringToMap(eqns);
        }
        return r1rhoEqnMap;
    }

    static Map<String, Boolean> getExpEqnMap() {
        if (expEqnMap == null) {
            Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
            String eqns = prefs.get("EXP_EQNS", DEFAULT_EXP_EQNS);
            expEqnMap = stringToMap(eqns);
        }
        return expEqnMap;
    }

    static Map<String, Boolean> getNOEEqnMap() {
        if (noeEqnMap == null) {
            Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
            String eqns = prefs.get("NOE_EQNS", DEFAULT_NOE_EQNS);
            noeEqnMap = stringToMap(eqns);
        }
        return noeEqnMap;
    }

    public static void setCESTEqnMap(String eqns) {
        cestEqnMap = stringToMap(eqns);
    }
    public static void setCPMGEqnMap(String eqns) {
        cpmgEqnMap = stringToMap(eqns);
    }

    public static void setR1RhoEqnMap(String eqns) {
        r1rhoEqnMap = stringToMap(eqns);
    }

    public static void setExpEqnMap(String eqns) {
        expEqnMap = stringToMap(eqns);
    }

    static Map<String, Boolean> stringToMap(String s) {
        Map<String, Boolean> map = new HashMap<>();
        String[] stringList = s.split("\n");
        for (String strValue : stringList) {
            String[] strValueParts = strValue.split(";");
            String key = strValueParts[0];
            if (strValueParts[1].equals("true")) {
                map.put(key, Boolean.TRUE);
            } else {
                map.put(key, Boolean.FALSE);
            }
        }
        return map;
    }

    static String mapToString(Map<String, Boolean> map) {
        StringBuilder sBuilder = new StringBuilder();
        for (String eqn : map.keySet()) {
            sBuilder.append(eqn);
            sBuilder.append(';');
            sBuilder.append(map.get(eqn));
            sBuilder.append("\n");
        }
        return sBuilder.toString();
    }

    public static void saveCESTEqnPrefs() {
        Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
        String eqnString = DEFAULT_CEST_EQNS;
        if (cestEqnMap != null) {
            eqnString = mapToString(getCESTEqnMap());
        }
        prefs.put("CEST_EQNS", eqnString);
    }
    public static void saveCPMGEqnPrefs() {
        Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
        String eqnString = DEFAULT_CPMG_EQNS;
        if (cpmgEqnMap != null) {
            eqnString = mapToString(getCPMGEqnMap());
        }
        prefs.put("CPMG_EQNS", eqnString);
    }

    public static void saveR1RhoEqnPrefs() {
        Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
        String eqnString = DEFAULT_R1RHO_EQNS;
        if (r1rhoEqnMap != null) {
            eqnString = mapToString(getR1RhoEqnMap());
        }
        prefs.put("R1RHO_EQNS", eqnString);
    }

    public static void saveExpEqnPrefs() {
        Preferences prefs = getPrefs();//Preferences.userNodeForPackage(ExperimentData.class);
        String eqnString = DEFAULT_EXP_EQNS;
        if (expEqnMap != null) {
            eqnString = mapToString(getExpEqnMap());
        }
        prefs.put("EXP_EQNS", eqnString);
    }

    public static List<String> getActiveCESTEquations() {
        Map<String, Boolean> map = getCESTEqnMap();
        List<String> cestEqnList = new ArrayList<>();
        for (String eqn : map.keySet()) {
            if (map.get(eqn)) {
                cestEqnList.add(eqn);
            }
        }
        return cestEqnList;
    }
    public static List<String> getActiveCPMGEquations() {
        Map<String, Boolean> map = getCPMGEqnMap();
        List<String> cpmgEqnList = new ArrayList<>();
        for (String eqn : map.keySet()) {
            if (map.get(eqn)) {
                cpmgEqnList.add(eqn);
            }
        }
        return cpmgEqnList;
    }

    public static List<String> getActiveR1RhoEquations() {
        Map<String, Boolean> map = getR1RhoEqnMap();
        List<String> r1rhoEqnList = new ArrayList<>();
        for (String eqn : map.keySet()) {
            if (map.get(eqn)) {
                r1rhoEqnList.add(eqn);
            }
        }
        return r1rhoEqnList;
    }

    public static List<String> getActiveExpEquations() {
        Map<String, Boolean> map = getExpEqnMap();
        List<String> expEqnList = new ArrayList<>();
        for (String eqn : map.keySet()) {
            if (map.get(eqn)) {
                expEqnList.add(eqn);
            }
        }
        return expEqnList;
    }

    public static List<String> getActiveNOEEquations() {
        Map<String, Boolean> map = getNOEEqnMap();
        List<String> noeEqnList = new ArrayList<>();
        for (String eqn : map.keySet()) {
            if (map.get(eqn)) {
                noeEqnList.add(eqn);
            }
        }
        return noeEqnList;
    }

    public static void setCESTEquationState(String equation, boolean state) {
        Map<String, Boolean> map = getCESTEqnMap();
        map.put(equation, state);
        saveCESTEqnPrefs();
    }

    public static boolean getCESTEquationState(String equation) {
        Map<String, Boolean> map = getCESTEqnMap();
        boolean state = false;
        if (map.containsKey(equation) && map.get(equation)) {
            state = true;
        }
        return state;
    }

    public static void setCPMGEquationState(String equation, boolean state) {
        Map<String, Boolean> map = getCPMGEqnMap();
        map.put(equation, state);
        saveCPMGEqnPrefs();
    }

    public static boolean getCPMGEquationState(String equation) {
        Map<String, Boolean> map = getCPMGEqnMap();
        boolean state = false;
        if (map.containsKey(equation) && map.get(equation)) {
            state = true;
        }
        return state;
    }

    public static void setR1RhoEquationState(String equation, boolean state) {
        Map<String, Boolean> map = getR1RhoEqnMap();
        map.put(equation, state);
        saveR1RhoEqnPrefs();
    }

    public static boolean getR1RhoEquationState(String equation) {
        Map<String, Boolean> map = getR1RhoEqnMap();
        boolean state = false;
        if (map.containsKey(equation) && map.get(equation)) {
            state = true;
        }
        return state;
    }

    public static void setExpEquationState(String equation, boolean state) {
        Map<String, Boolean> map = getExpEqnMap();
        map.put(equation, state);
        saveExpEqnPrefs();
    }

    public static boolean getExpEquationState(String equation) {
        Map<String, Boolean> map = getExpEqnMap();
        boolean state = false;
        if (map.containsKey(equation) && map.get(equation)) {
            state = true;
        }
        return state;
    }

}
