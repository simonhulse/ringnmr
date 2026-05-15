package org.comdnmr.modelfree;

import org.comdnmr.data.DynamicsSource;
import org.nmrfx.chemistry.Atom;
import java.util.ArrayList;
import java.util.List;

public class R1R2NOEMolDataValues extends MolDataValues<R1R2NOEDataValue> {

    /**
     * Strategy for computing J(0) from multi-field relaxation data.
     *
     * <ul>
     *   <li>{@link #INDEPENDENT} — one J(0) per field, giving a 3F J-value
     *       vector (F = number of fields). This is the default and matches
     *       the original implementation.</li>
     *   <li>{@link #AVERAGED_JACKKNIFE} — a single J(0) derived from a
     *       no-intercept weighted least-squares regression of the Γ values
     *       across all fields, with jackknife uncertainty. The resulting
     *       J-value vector has length 1 + 2F.</li>
     * </ul>
     */
    public enum J0Mode {
        INDEPENDENT,
        AVERAGED_JACKKNIFE
    }

    private J0Mode j0Mode = J0Mode.INDEPENDENT;

    /** Per-field Γ regression weights; only used in {@link J0Mode#AVERAGED_JACKKNIFE} mode. */
    private double[] j0Weights = null;

    public R1R2NOEMolDataValues(String specifier, double[] vector, DynamicsSource dynSourceFactory) {
        super(specifier, vector, dynSourceFactory);
    }
    public R1R2NOEMolDataValues(Atom atom, double[] vector) { super(atom, vector); }
    public R1R2NOEMolDataValues(Atom atom)                  { super(atom); }

    public J0Mode getJ0Mode() { return j0Mode; }

    public void setJ0Mode(J0Mode mode) {
        this.j0Mode = mode;
        clearJValues();
    }

    public void setJ0Weights(double[] weights) {
        this.j0Weights = weights;
        clearJValues();
    }

    @Override
    public double[][] calcJ() {
        if (dataValues.isEmpty()) return new double[0][0];
        return j0Mode == J0Mode.AVERAGED_JACKKNIFE
            ? SpectralDensityCalculator.calcJR1R2NOEAveraged(dataValues, j0Weights)
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
        return j0Mode == J0Mode.AVERAGED_JACKKNIFE
            ? 1 + dataValues.size() * 2
            : dataValues.size() * 3;
    }

    @Override
    public MolDataValues<R1R2NOEDataValue> createEmpty() {
        return new R1R2NOEMolDataValues(atom, vector);
    }

    @Override
    public R1R2NOEMolDataValues copy() {
        R1R2NOEMolDataValues copy = new R1R2NOEMolDataValues(atom, vector);
        copy.setTestModel(getTestModel());
        copy.j0Mode = this.j0Mode;
        if (j0Weights != null) copy.j0Weights = j0Weights.clone();
        for (R1R2NOEDataValue dv : dataValues) {
            copy.addData(new R1R2NOEDataValue(copy, dv.R1, dv.R1err, dv.R2, dv.R2err, dv.NOE, dv.NOEerr, dv.relaxObj));
        }
        return copy;
    }
}
