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
    public AmideBayesianSampler(MolDataValues<R1R2NOEDataValue> data) {
        super(data);
        dirichlets = new DirichletSampler[N_FREQ_CLASSES];
        for (int i = 0; i < N_FREQ_CLASSES; i++) {
            dirichlets[i] = DirichletSampler.symmetric(rng, getNFields(), ALPHA);
        }
    }

    /**
     * Draws per-frequency-class Dirichlet weights and assembles them into a flat weight
     * vector using the same index layout as {@link AmideNonparametricSampler}:
     * {@code weights[fieldIndex * 3 + freqClass]}.
     *
     * @return per-observation bootstrap weights
     */
    @Override
    protected double[] sampleWeights() {
        int nFields = getNFields();
        double[] weights = new double[getNSpectralDensities()];
        for (int i = 0; i < N_FREQ_CLASSES; i++) {
            double[] w = dirichlets[i].sample();
            for (int k = 0; k < nFields; k++) {
                weights[k * 3 + i] = w[k] * nFields;
            }
        }
        return weights;
    }
}
