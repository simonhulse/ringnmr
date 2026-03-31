package org.comdnmr.modelfree;

import org.apache.commons.rng.sampling.distribution.DirichletSampler;

public class BayesianSampler extends WeightSampler {

    // Concentration factor
    private static final double ALPHA = 1.0;

    // Dirichlet sampler instance
    private final DirichletSampler dirichlet;

    public BayesianSampler(MolDataValues data) {
        super(data);
        dirichlet = DirichletSampler.symmetric(rng, getNValues(), ALPHA);
    }

    public double[] sampleWeights() {
        int nJ = getNValues();
        double[] weights = dirichlet.sample();
        // Weights for Dirichlet distribution sum to 1.
        // Scale these so that they sum to nJ.
        for (int j = 0; j < nJ; j++) weights[j] *= nJ;
        return weights;
    }
}
