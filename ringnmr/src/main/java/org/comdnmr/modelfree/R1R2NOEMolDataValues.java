package org.comdnmr.modelfree;

import org.comdnmr.data.DynamicsSource;
import org.nmrfx.chemistry.Atom;
import java.util.ArrayList;
import java.util.List;

public class R1R2NOEMolDataValues extends MolDataValues<R1R2NOEDataValue> {

    public R1R2NOEMolDataValues(String specifier, double[] vector, DynamicsSource dynSourceFactory) {
        super(specifier, vector, dynSourceFactory);
    }
    public R1R2NOEMolDataValues(Atom atom, double[] vector) { super(atom, vector); }
    public R1R2NOEMolDataValues(Atom atom)                  { super(atom); }

    @Override
    public double[][] calcJ() {
        return dataValues.isEmpty() ? new double[0][0]
            : SpectralDensityCalculator.calcJR1R2NOE(dataValues);
    }

    @Override
    public List<double[][]> calcIndependentJ() {
        List<double[][]> result = new ArrayList<>();
        for (R1R2NOEDataValue value : dataValues) {
            result.add(SpectralDensityCalculator.calcJR1R2NOE(List.of(value)));
        }
        return result;
    }

    @Override
    public int getNSpectralDensities() {
        return dataValues.size() * 3;
    }

    @Override
    public MolDataValues<R1R2NOEDataValue> createEmpty() {
        return new R1R2NOEMolDataValues(atom, vector);
    }

    @Override
    public R1R2NOEMolDataValues copy() {
        R1R2NOEMolDataValues copy = new R1R2NOEMolDataValues(atom, vector);
        copy.setTestModel(getTestModel());
        for (R1R2NOEDataValue dv : dataValues) {
            copy.addData(new R1R2NOEDataValue(copy, dv.R1, dv.R1err, dv.R2, dv.R2err, dv.NOE, dv.NOEerr, dv.relaxObj));
        }
        return copy;
    }
}
