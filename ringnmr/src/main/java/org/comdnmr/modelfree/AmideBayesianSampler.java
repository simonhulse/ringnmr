package org.comdnmr.modelfree;

import org.apache.commons.rng.sampling.distribution.DirichletSampler;

/**
 * Bayesian bootstrap sampler that applies the same segmented weighting strategy as
 * {@link AmideNonparametricSampler}, but draws continuous weights from independent
 * Dirichlet distributions rather than discrete selection counts.
 *
 * <p>For the three canonical amide spectral-density frequency classes — J(0), J(ω_N),
 * J(0.87ω_H) — a separate Dirichlet(1, …, 1) sampler of dimension {@code nFields} is
 * used. Each draw is scaled so that the weights for a given frequency class sum to
 * {@code nFields}, matching the scale of the integer counts produced by
 * {@link AmideNonparametricSampler}.
 *
 * <p>Because the Dirichlet distribution is continuous, this sampler can be called an
 * unlimited number of times, unlike {@link AmideNonparametricSampler}.
 */
public class AmideBayesianSampler extends WeightSampler<R1R2NOEDataValue> {

    private static final double ALPHA = 1.0;
    private static final int N_FREQ_CLASSES = 3;

    private final DirichletSampler[] dirichlets;

    /**
     * Constructs an {@code AmideBayesianSampler} for the given relaxation data.
     *
     * @param data the relaxation data to be resampled
     */
    public AmideBayesianSampler(MolDataValues<R1R2NOEDataValue> data) { this(data, false); }

    // Added for use in the regularization paper; not used within RING.
    public static AmideBayesianSampler withFixedSeed(MolDataValues<R1R2NOEDataValue> data) {
        return new AmideBayesianSampler(data, true);
    }

    private AmideBayesianSampler(MolDataValues<R1R2NOEDataValue> data, boolean seed) {
        super(data, seed);
        dirichlets = new DirichletSampler[N_FREQ_CLASSES];
        for (int i = 0; i < N_FREQ_CLASSES; i++) {
            dirichlets[i] = DirichletSampler.symmetric(rng, getNFields(), ALPHA);
        }
    }

    /**
     * Returns the next bootstrap replicate, branching on the data's
     * {@link R1R2NOEMolDataValues.J0Mode}.
     *
     * <p>For {@link R1R2NOEMolDataValues.J0Mode#INDEPENDENT} the copy carries a
     * 3F weight vector (Dirichlet draws scaled by F per frequency class).
     *
     * <p>For {@link R1R2NOEMolDataValues.J0Mode#AVERAGED_JACKKNIFE} the copy
     * carries per-field Γ regression weights (via
     * {@link R1R2NOEMolDataValues#setJ0Weights}) and a (1 + 2F) chi-sq weight
     * vector with 1.0 at index 0 for the single J(0).
     *
     * @return a new {@link MolDataValues} representing the current bootstrap replicate
     */
    @Override
    public MolDataValues<R1R2NOEDataValue> sample() {
        int nFields = getNFields();
        double[] d0 = dirichlets[0].sample();
        double[] d1 = dirichlets[1].sample();
        double[] d2 = dirichlets[2].sample();

        R1R2NOEMolDataValues copy = (R1R2NOEMolDataValues) data.copy();

        if (((R1R2NOEMolDataValues) data).getJ0Mode() == R1R2NOEMolDataValues.J0Mode.AVERAGED_JACKKNIFE) {
            double[] j0Weights = new double[nFields];
            for (int k = 0; k < nFields; k++) j0Weights[k] = d0[k] * nFields;
            copy.setJ0Weights(j0Weights);

            double[] chiSqWeights = new double[1 + 2 * nFields];
            chiSqWeights[0] = 1.0;
            for (int k = 0; k < nFields; k++) {
                chiSqWeights[1 + k * 2]     = d1[k] * nFields;
                chiSqWeights[1 + k * 2 + 1] = d2[k] * nFields;
            }
            copy.setWeights(chiSqWeights);
        } else {
            double[] weights = new double[3 * nFields];
            for (int k = 0; k < nFields; k++) {
                weights[k * 3]     = d0[k] * nFields;
                weights[k * 3 + 1] = d1[k] * nFields;
                weights[k * 3 + 2] = d2[k] * nFields;
            }
            copy.setWeights(weights);
        }
        return copy;
    }
}
