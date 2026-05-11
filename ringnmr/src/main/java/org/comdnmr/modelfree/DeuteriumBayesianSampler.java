package org.comdnmr.modelfree;

import org.apache.commons.rng.sampling.distribution.DirichletSampler;

/**
 * Bayesian bootstrap sampler for deuterium relaxation data.
 *
 * <p>Applies the same segmented weighting strategy as {@link DeuteriumNonparametricSampler}
 * but draws continuous Dirichlet weights rather than discrete counts. J(0) is always
 * assigned a constant weight of 1. The J(ω_D) and J(2ω_D) groups each receive an
 * independent Dirichlet(1, …, 1) draw, scaled to sum to the group size.
 *
 * <p>Degenerate frequencies (where ω_D of one field coincides with 2ω_D of another)
 * are placed in the J(ω_D) group; see {@link DeuteriumFrequencyGroups}.
 *
 * <p>This sampler can be called an unlimited number of times, unlike
 * {@link DeuteriumNonparametricSampler}.
 */
public class DeuteriumBayesianSampler extends WeightSampler<DeuteriumDataValue> {

    private static final double ALPHA = 1.0;

    private final DeuteriumFrequencyGroups groups;
    private final DirichletSampler wdDirichlet;
    private final DirichletSampler w2dDirichlet;

    /**
     * Constructs a {@code DeuteriumBayesianSampler} for the given relaxation data.
     *
     * @param data the deuterium relaxation data to be resampled
     */
    public DeuteriumBayesianSampler(MolDataValues<DeuteriumDataValue> data) {
        super(data);
        groups       = DeuteriumFrequencyGroups.of(data.getData());
        wdDirichlet  = groups.wdIndices.length  > 0 ? DirichletSampler.symmetric(rng, groups.wdIndices.length,  ALPHA) : null;
        w2dDirichlet = groups.w2dIndices.length > 0 ? DirichletSampler.symmetric(rng, groups.w2dIndices.length, ALPHA) : null;
    }

    /**
     * Draws per-group Dirichlet weights. J(0) is fixed at 1; J(ω_D) and J(2ω_D)
     * observations receive weights drawn from Dirichlet(1, …, 1), scaled to sum to
     * their respective group sizes.
     *
     * @return per-observation bootstrap weights
     */
    @Override
    protected double[] sampleWeights() {
        double[] weights = new double[getNSpectralDensities()];
        weights[0] = 1.0;
        if (wdDirichlet != null) {
            double[] w = wdDirichlet.sample();
            int n = groups.wdIndices.length;
            for (int i = 0; i < n; i++) weights[groups.wdIndices[i]] = w[i] * n;
        }
        if (w2dDirichlet != null) {
            double[] w = w2dDirichlet.sample();
            int n = groups.w2dIndices.length;
            for (int i = 0; i < n; i++) weights[groups.w2dIndices[i]] = w[i] * n;
        }
        return weights;
    }
}
