package org.comdnmr.modelfree;

import java.util.Arrays;

/**
 * Abstract base for weight-based bootstrap samplers.
 *
 * <p>Instead of perturbing observed values, a {@code WeightSampler} assigns a
 * non-negative weight to each spectral-density observation. Concrete subclasses
 * implement {@link #sample()} to produce a copy of the wrapped data with
 * appropriate weights set via {@link MolDataValues#setWeights(double[])}.
 *
 * @see AmideNonparametricSampler
 * @see AmideBayesianSampler
 * @see DeuteriumNonparametricSampler
 * @see DeuteriumBayesianSampler
 */
public abstract class WeightSampler<T extends RelaxDataValue> extends BootstrapSampler<T> {

    WeightSampler(MolDataValues<T> data) { super(data); }

    WeightSampler(MolDataValues<T> data, boolean seed) { super(data, seed); }

    /**
     * {@inheritDoc}
     */
    public MolDataValues<T> getOriginalData() {
        return data;
    }

}
