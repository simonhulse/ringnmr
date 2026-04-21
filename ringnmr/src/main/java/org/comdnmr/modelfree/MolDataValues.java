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
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

/**
 *
 * @author brucejohnson
 */
public abstract class MolDataValues<T extends RelaxDataValue> {

    Atom atom;
    final String specifier;
    final double[] vector = new double[3];
    final List<T> dataValues = new ArrayList<>();
    private MFModel model;
    private double[][] jValues = null;
    private Integer bootstrapSet = null;
    private BootstrapAggregator bootstrapAggregator = null;
    private double[] weights = null;

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

    public String getSpecifier() {
        return specifier;
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

    public void setWeights(double[] weights) {
        this.weights = weights;
    }

    /** @deprecated Use {@link #setWeights(double[])} instead. */
    @Deprecated
    public void weight(double[] weights) {
        setWeights(weights);
    }

    public double[] getWeights() {
        double[] weights;
        if (this.weights == null) {
            int nWeights = getNValues();
            weights = new double[nWeights];
            Arrays.fill(weights, 1.0);
        } else {
            weights = this.weights;
        }
        return weights;
    }

    public void addData(T value) {
        dataValues.add(value);
    }

    public List<T> getData() {
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

    public void clearJValues() {
        jValues = null;
    }

    public abstract double[][] calcJ();

    public abstract List<double[][]> calcIndependentJ();

    public abstract MolDataValues<T> createEmpty();

    public void setJValues(double[][] jValuesSet) {
        jValues = Arrays.stream(jValuesSet).map(double[]::clone).toArray(double[][]::new);
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
            System.arraycopy(weights, 0, jValues[jValues.length - 1], 0, weights.length);
        }
        return jValues;
    }
}
