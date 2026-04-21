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

    @Override
    public MolDataValues<DeuteriumDataValue> createEmpty() {
        return new DeuteriumMolDataValues(atom, vector);
    }
}
