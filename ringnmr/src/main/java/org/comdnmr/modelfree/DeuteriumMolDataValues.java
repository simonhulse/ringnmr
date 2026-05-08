package org.comdnmr.modelfree;

import org.comdnmr.data.DynamicsSource;
import org.nmrfx.chemistry.Atom;
import java.util.ArrayList;
import java.util.List;

public class DeuteriumMolDataValues extends MolDataValues<DeuteriumDataValue> {

    public DeuteriumMolDataValues(String specifier, double[] vector, DynamicsSource dynSourceFactory) {
        super(specifier, vector, dynSourceFactory);
    }
    public DeuteriumMolDataValues(Atom atom, double[] vector) { super(atom, vector); }
    public DeuteriumMolDataValues(Atom atom)                  { super(atom); }

    @Override
    public double[][] calcJ() {
        return dataValues.isEmpty() ? new double[0][0]
            : SpectralDensityCalculator.calcJDeuterium(dataValues);
    }

    @Override
    public List<double[][]> calcIndependentJ() {
        List<double[][]> result = new ArrayList<>();
        for (DeuteriumDataValue value : dataValues) {
            result.add(SpectralDensityCalculator.calcJDeuterium(List.of(value)));
        }
        return result;
    }

    private Integer cachedNSpectralDensities = null;

    @Override
    public void addData(DeuteriumDataValue value) {
        cachedNSpectralDensities = null;
        super.addData(value);
    }

    @Override
    public int getNSpectralDensities() {
        if (cachedNSpectralDensities != null) return cachedNSpectralDensities;
        if (dataValues.isEmpty()) return 0;
        List<Double> fieldList = new ArrayList<>();
        fieldList.add(0.0);
        for (DeuteriumDataValue dv : dataValues) {
            double omega = dv.getB0() * RelaxEquations.GAMMA_D / RelaxEquations.GAMMA_H * 2.0 * Math.PI;
            if (fieldList.stream().noneMatch(f -> f > 0 && Math.abs((omega - f) / f) < 0.01)) {
                fieldList.add(omega);
            }
            double omega2 = omega * 2.0;
            if (fieldList.stream().noneMatch(f -> f > 0 && Math.abs((omega2 - f) / f) < 0.01)) {
                fieldList.add(omega2);
            }
        }
        cachedNSpectralDensities = fieldList.size();
        return cachedNSpectralDensities;
    }

    @Override
    public MolDataValues<DeuteriumDataValue> createEmpty() {
        return new DeuteriumMolDataValues(atom, vector);
    }
}
