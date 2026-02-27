/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.gui;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.comdnmr.data.ExperimentResult;
import org.comdnmr.data.ExperimentSet;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.ValueSet;

/**
 *
 * @author brucejohnson
 */
public class ChartInfo {

    List<String> mapName = new ArrayList<>();
    String state;
    String equationName;
    ResidueChart.SelectionValue[] currentSelections;
    ExperimentResult experimentalResult;
    ValueSet valueSet;
    List<int[]> currentStates = new ArrayList<>();

    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(mapName).append(" ").append(equationName).append(" ").append(state).append(" ").
                append(valueSet.name());
        for (var selectionValue : currentSelections) {
            sBuilder.append(" ").append(selectionValue.resonanceSource().getAtom().getShortName());
        }
        return sBuilder.toString();
    }

    public boolean hasExperiments() {
        return valueSet != null && (valueSet instanceof ExperimentSet);
    }

    public ValueSet getExperiments() {
        return valueSet;
    }

    public boolean hasResult() {
        return experimentalResult != null;
    }

    public boolean hasResidue() {
        return experimentalResult != null;
    }

    public boolean hasResidues() {
        return currentSelections != null;
    }

    public Atom getAtom() {
        return experimentalResult.getAtom();
    }

    public int getResNum() {
        return experimentalResult.getAtom().getResidueNumber();
    }

    public ResonanceSource getSource() {
        return experimentalResult.getResonanceSource();
    }

    public ResidueChart.SelectionValue[] getCurrentSelections() {
        return currentSelections;
    }

    public List<ResonanceSource> getCurrentResonanceSources() {
        return  Arrays.stream(currentSelections).map(sV -> sV.resonanceSource()).distinct().toList();
    }

    public void setResidues(ResidueChart.SelectionValue[] residues) {
        this.currentSelections = residues;
    }

    public ExperimentResult getResult() {
        return experimentalResult;
    }

    public String getEquationName() {
        return equationName;
    }

    public void setStates(List<int[]> states) {
        currentStates.clear();
        currentStates.addAll(states);
    }

    public void clear() {
        currentSelections = null;
        experimentalResult = null;
        valueSet = null;
        currentStates = new ArrayList<>();
    }

}
