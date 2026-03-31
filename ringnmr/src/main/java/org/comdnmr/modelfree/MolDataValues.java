/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import org.comdnmr.data.DynamicsSource;
import org.comdnmr.modelfree.models.MFModel;
import org.nmrfx.chemistry.Atom;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

/**
 *
 * @author brucejohnson
 */
public class MolDataValues {

    Atom atom;
    final String specifier;
    final double[] vector = new double[3];
    final List<RelaxDataValue> dataValues = new ArrayList<>();
    MFModel model;
    double[][] jValues = null;
    Integer bootstrapSet = null;
    BootstrapAggregator bootstrapAggregator = null;
    double[] weights = null;

    public MolDataValues(String specifier, double[] vector, DynamicsSource dynSourceFactory) {
        this.specifier = specifier;
        Optional<Atom> atomOpt = dynSourceFactory.getAtom(null, specifier, 0, false);
        atom = atomOpt.orElse(null);
        System.arraycopy(vector, 0, this.vector, 0, 3);
    }

    public MolDataValues(Atom atom, double[] vector) {
        this.atom = atom;
        this.specifier = atom.getFullName();
        System.arraycopy(vector, 0, this.vector, 0, 3);
    }

    public MolDataValues(Atom atom) {
        this.atom = atom;
        this.specifier = atom.getFullName();
    }

    public Atom getAtom() {
        return atom;
    }
    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        for (var dataValue : dataValues) {
            sBuilder.append(dataValue.toString());
            sBuilder.append("\n");
        }
        return sBuilder.toString();
    }

    public int getNValues() {
        return SpectralDensityCalculator.getNData(dataValues);
    }

    public void weight(double[] weights) {
        this.weights = weights;
    }

    public double[] getWeights() {
        double[] weights;
        if (this.weights == null) {
            int nWeights = getNValues();
            weights = new double[nWeights];
            for (int j = 0; j < nWeights; j++) {
                weights[j] = 1.0;
            }
        } else {
            weights = this.weights;
        }
        return weights;
    }

    public void addData(RelaxDataValue value) {
        dataValues.add(value);
    }

    public List<RelaxDataValue> getData() {
        return dataValues;
    }

    public void setTestModel(MFModel model) {
        this.model = model;
    }

    public MFModel getTestModel() {
        return model;
    }

    public double[] getVector() {
        return vector;
    }

    public void setBootstrapAggregator(BootstrapAggregator bootstrapAggregator) {
        this.bootstrapAggregator = bootstrapAggregator;
    }

    public void setBootstrapSet(int iSet) {
        bootstrapSet = iSet;
        jValues = null;
    }

    public void clearBootStrapSet() {
        bootstrapSet = null;
        jValues = null;
    }

    public void clearJ() {
        jValues = null;
    }

    public double[][] calcJ() {
        var dataOpt = dataValues.stream().findFirst();
        if (dataOpt.isPresent()) {
            if (dataOpt.get() instanceof DeuteriumDataValue) {
                return SpectralDensityCalculator.calcJDeuterium(dataValues);
            } else {
                return SpectralDensityCalculator.calcJR1R2NOE(dataValues);
            }
        }
        return new double[0][0];
    }

    public List<double[][]> calcIndependentJ() {
        var dataOpt = dataValues.stream().findFirst();
        List<double[][]> result = new ArrayList<>();
        if (dataOpt.isPresent()) {
            if (dataOpt.get() instanceof DeuteriumDataValue) {
                for (RelaxDataValue relaxDataValue : dataValues) {
                    double[][] jValues = SpectralDensityCalculator.calcJDeuterium(List.of(relaxDataValue));
                    result.add(jValues);
                }
            } else {
                for (RelaxDataValue relaxDataValue : dataValues) {
                    double[][] jValues = SpectralDensityCalculator.calcJR1R2NOE(List.of(relaxDataValue));
                    result.add(jValues);
                }
                result.add(jValues);
            }
        }
        return result;
    }

    public void setJValues(double[][] jValuesSet) {
        jValues = new double[jValuesSet.length][jValuesSet[0].length];
        for (int i=0;i<jValuesSet.length;i++) {
            for (int j=0;j<jValuesSet[i].length;j++) {
                jValues[i][j] = jValuesSet[i][j];
            }
        }
    }

    public double[][] getJValues() {
        if (jValues == null) {
            if (bootstrapSet != null) {
                jValues = bootstrapAggregator.getBootStrapJ(calcJ(), bootstrapSet);
            } else {
                jValues = calcJ();
            }
        }
        if (weights != null) {
            System.arraycopy(weights, 0, jValues[jValues.length -1], 0, weights.length);
        }
        return jValues;

    }
}
