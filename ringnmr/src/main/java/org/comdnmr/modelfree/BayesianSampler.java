package org.comdnmr.modelfree;

import org.apache.commons.rng.sampling.distribution.DirichletSampler;

/**
 * Bayesian bootstrap sampler that draws observation weights from a symmetric
 * Dirichlet distribution.
 *
 * <p>The Bayesian bootstrap (Rubin 1981) replaces the classical bootstrap's
 * integer-valued resampling counts with continuous weights drawn from
 * Dirichlet(α, …, α) with α = 1 (the uniform Dirichlet prior). Each call to
 * {@link #sampleWeights()} returns a weight vector whose elements sum to
 * {@code nValues} rather than 1, so that the downstream weighted least-squares
 * objective is on the same scale as the unweighted case.
 *
 * <p>This sampler can be called an unlimited number of times, unlike
 * {@link NonparametricSampler} which exhausts after {@code nSel³} draws.
 */
public class BayesianSampler<T extends RelaxDataValue> extends WeightSampler<T> {

    /** Dirichlet concentration parameter; α = 1 corresponds to the uniform prior. */
    private static final double ALPHA = 1.0;

    private final DirichletSampler dirichlet;

    /**
     * Constructs a {@code BayesianSampler} for the given relaxation data.
     *
     * @param data the relaxation data to be resampled
     */
    public BayesianSampler(MolDataValues<T> data) {
        super(data);
        dirichlet = DirichletSampler.symmetric(rng, getNSpectralDensities(), ALPHA);
    }

    /**
     * Draws a weight vector from Dirichlet(1, …, 1) and scales it so that the
     * weights sum to {@code nValues}.
     *
     * @return per-observation bootstrap weights summing to {@code getNSpectralDensities()}
     */
    public double[] sampleWeights() {
        int nJ = getNSpectralDensities();
        double[] weights = dirichlet.sample();
        // Weights for Dirichlet distribution sum to 1.
        // Scale these so that they sum to nJ.
        for (int j = 0; j < nJ; j++) weights[j] *= nJ;
        return weights;
    }
}
