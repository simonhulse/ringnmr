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
package org.comdnmr.data;

import java.util.*;
import java.util.stream.Collectors;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.relax.RelaxationValues;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.ValueSet;
import org.nmrfx.chemistry.relax.ValueWithError;

/**
 *
 * @author Bruce Johnson
 */
public class ExperimentSet implements ValueSet {

    Map<String, Experiment> expMaps = new HashMap<>();
    private final HashMap<String, ExperimentResult> resultMap = new HashMap<>();
    final String name;
    String fileName = null;
    Map<Double, Integer> fieldMap = new LinkedHashMap();
    Map<Double, Integer> tempMap = new LinkedHashMap();
    Map<Double, Integer> b1Map = new LinkedHashMap();
    List<Double> fieldList = new ArrayList<>();
    List<Double> tempList = new ArrayList<>();
    List<Double> b1List = new ArrayList<>();
    private boolean absValueMode = false;
    private String bootStrapMode = "parametric";
    private String expMode = "cpmg";

    private boolean active = true;

    public ExperimentSet(String name, String fileName) {
        this.name = name;
        this.fileName = fileName;
    }

    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(name).append(" ").append(fileName).append(" resmap ").append(resultMap.size());
        return sBuilder.toString();
    }

    public String name() {
        return name;
    }

    /**
     * @return the absValueMode
     */
    public boolean isAbsValueMode() {
        return absValueMode;
    }

    /**
     * @param absValueMode the absValueMode to set
     */
    public void setAbsValueMode(boolean absValueMode) {
        this.absValueMode = absValueMode;
    }

    /**
     * @return the bootStrapMode
     */
    public String getBootStrapMode() {
        return bootStrapMode;
    }

    /**
     * @param bootStrapMode the bootStrapMode to set
     */
    public void setBootStrapMode(String bootStrapMode) {
        this.bootStrapMode = bootStrapMode;
    }

    /**
     * @return the expMode
     */
    public String getExpMode() {
        return expMode;
    }

    /**
     * @param expMode the expMode to set
     */
    public void setExpMode(String expMode) {
        this.expMode = expMode;
    }

    public void addExperimentData(String name, Experiment data) {
        expMaps.put(name, data);
    }

    public Experiment getExperimentData(String name) {
        return expMaps.get(name);
    }

    public Collection<Experiment> getExperimentData() {
        return expMaps.values();
    }

    public Map<String, Experiment> getExperimentMap() {
        return expMaps;
    }

    public List<String> getEquationNames() {
        Set<String> equationSet = new HashSet<>();
        resultMap.values().stream().forEach(rInfo -> {
            equationSet.addAll(rInfo.curveFits.keySet());
        });
        List<String> equationNames = equationSet.stream().sorted().collect(Collectors.toList());
        return equationNames;
    }

    public static boolean matchStateString(String target, String state) {
        boolean match = true;
        String[] targetElems = target.split(":");
        String[] stateElems = state.split(":");
        for (int i = 0; i < targetElems.length; i++) {
            if (targetElems.length == stateElems.length) {
                if (!targetElems[i].equals("*") && !targetElems[i].equals(stateElems[i])) {
                    match = false;
                    break;
                }
            }
        }
        return match;
    }

    public List<String> getStateStrings() {
        if (fieldMap.isEmpty()) {
            setupMaps();
        }
        // List<String> stateStrings = 
        Set<String> states = expMaps.values().stream().map(e -> e.getState()).map(s
                -> "0:" + s.substring(2)).collect(Collectors.toSet());
        List<String> stateStrings = states.stream().sorted().collect(Collectors.toList());
        return stateStrings;
    }

    public static String getStateString(int[] state) {
        StringBuilder builder = new StringBuilder();
        builder.append(state[1]);
        for (int i = 2; i < state.length; i++) {
            builder.append(':');
            builder.append(state[i]);
        }
        return builder.toString();
    }

    public void setupMaps() {
        fieldMap.clear();
        tempMap.clear();
        b1Map.clear();
        fieldList.clear();
        tempList.clear();
        b1List.clear();
        for (Experiment experiment : expMaps.values()) {
            if (!fieldMap.containsKey(Math.floor(experiment.getB0Field()))) {
                fieldMap.put(Math.floor(experiment.getB0Field()), fieldMap.size());
                fieldList.add(experiment.getB0Field());
            }
            if (!tempMap.containsKey(Math.floor(experiment.getTemperature()))) {
                tempMap.put(Math.floor(experiment.getTemperature()), tempMap.size());
                tempList.add(experiment.getTemperature());
            }
            double b1Field = 0.0;
            if (experiment instanceof OffsetExperiment offsetExperiment) {
                b1Field = offsetExperiment.getB1Field();
            }
            double b1Floor = Math.floor(b1Field);
            if (!b1Map.containsKey(b1Floor)) {
                b1Map.put(b1Floor, b1Map.size());
                b1List.add(b1Field);
            }
        }
        for (Experiment expData : expMaps.values()) {
            int[] state = getStateIndices(0, expData);
            expData.setState(getStateString(state));
        }
    }

    public int[] getStateIndices(int resIndex, Experiment experiment) {
        if (fieldMap.isEmpty()) {
            setupMaps();
        }
        int[] state = new int[4];
        state[0] = resIndex;
        state[1] = fieldMap.get(Math.floor(experiment.getB0Field()));
        state[2] = tempMap.get(Math.floor(experiment.getTemperature()));
        double b1Field = 0.0;
        if (experiment instanceof OffsetExperiment offsetExperiment) {
            b1Field = offsetExperiment.getB1Field();
        }
        double b1Floor = Math.floor(b1Field);

        state[3] = b1Map.get(b1Floor);
//        System.out.println(resIndex + " " + experiment.field + " " + experiment.temperature + " " + experiment.nucleus);
//        System.out.println("state index residue:" + state[0] + " field:" + state[1] + " temp:" + state[2] + " nuc:" + state[3]);

        return state;
    }

    public int[] getStateCount(int nResidues) {
        if (fieldMap.isEmpty()) {
            setupMaps();
        }
        int[] state = new int[4];
        state[0] = nResidues;
        state[1] = fieldMap.size();
        state[2] = tempMap.size();
        state[3] = b1Map.size();
//        System.out.println("state count residues:" + state[0] + " fields:" + state[1] + " temps:" + state[2] + " nucs:" + state[3]);
        return state;
    }

    public void addExperimentResult(ResonanceSource dynSource, ExperimentResult value) {
        resultMap.put(dynSource.getAtomKey(), value);
    }

    public void clearResidueMap() {
        resultMap.clear();
    }

    public Set<ResonanceSource> resonanceSources() {
        return resultMap.values().stream().map(v -> v.getResonanceSource()).filter(resSource -> !resSource.deleted()).collect(Collectors.toSet());
    }

    @Override
    public List<? extends RelaxationValues> rValues() {
        return Collections.emptyList();
    }

    @Override
    public void active(boolean b) {
        active = b;
    }

    @Override
    public boolean active() {
        return active;
    }

    public List<ExperimentResult> getExperimentResults() {
        List<ExperimentResult> values = new ArrayList<>();
        values.addAll(resultMap.values());
        return values;
    }

    public ExperimentResult getExperimentResult(ResonanceSource dynSource) {
        return resultMap.get(dynSource.getAtomKey());
    }

    public int getDataCount(ResonanceSource[] dynSources) {
        int n = 0;
        for (var dynSource : dynSources) {
            for (Experiment expData : expMaps.values()) {
                ExperimentData experimentalData = expData.getResidueData(dynSource);
                if (experimentalData != null) {
                    n++;
                }
            }
        }
        return n;
    }

    public double[] getFields() {
        double[] fields = new double[fieldList.size()];
        int i = 0;
        for (Double field : fieldList) {
            fields[i++] = field;
        }
        return fields;
    }

    public double[] getTemperatures() {
        double[] temperatures = new double[tempList.size()];
        int i = 0;
        for (Double temperature : tempList) {
            temperatures[i++] = temperature;
        }
        return temperatures;
    }

    public double[] getB1Field() {
        double[] b1Fields = new double[b1List.size()];
        int i = 0;
        for (double b1Field : b1List) {
            b1Fields[i++] = b1Field;
        }
        return b1Fields;
    }

    public Map<Atom, ValueWithError> getParMapData(String eqnName, String state, String parName) {

        List<ExperimentResult> resValues = getExperimentResults();
        Map<Atom, ValueWithError> resultMap = new HashMap<>();

        for (ExperimentResult experimentResult : resValues) {
            if (experimentResult == null) {
                continue;
            }
            if (experimentResult.getResonanceSource().deleted()) {
                continue;
            }
            String useEquName = eqnName;
            if (eqnName.equals("best")) {
                useEquName = experimentResult.getBestEquationName();
            }
            Atom atom = experimentResult.getAtom();
            Double y = experimentResult.getParValue(useEquName, state, parName);
            if (y == null) {
                continue;
            }
            Double yErr = experimentResult.getParValue(useEquName, state, parName+".sd");
            ValueWithError valueError = new ValueWithError(y, yErr);
            resultMap.put(atom, valueError);
        }
        return resultMap;
    }

}
