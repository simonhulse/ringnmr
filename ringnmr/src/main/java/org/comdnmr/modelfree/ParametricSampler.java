package org.comdnmr.modelfree;

import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratSampler;

/**
 * Parametric bootstrap sampler that generates synthetic observations by adding
 * Gaussian noise to the original measured values.
 *
 * <p>On each call to {@link #sample()}, a copy of the wrapped data is created and
 * every relaxation observable (R1, R2, and either NOE for {@link R1R2NOEDataValue}
 * or rQ/rAP for {@link DeuteriumDataValue}) is independently perturbed by a
 * zero-mean Gaussian with standard deviation equal to the corresponding reported
 * experimental error. The original data is never modified.
 *
 * <p>Normal deviates are drawn using the Ziggurat algorithm via Apache Commons RNG.
 */
public class ParametricSampler<T extends RelaxDataValue> extends BootstrapSampler<T> {

    private final NormalizedGaussianSampler gaussian;

    /**
     * Constructs a {@code ParametricSampler} for the given relaxation data.
     *
     * @param data the relaxation data to be resampled; must contain at least one entry
     */
    public ParametricSampler(MolDataValues<T> data) { this(data, false); }

    // Added for use in the regularization paper; not used within RING.
    public static <T extends RelaxDataValue> ParametricSampler<T> withFixedSeed(MolDataValues<T> data) {
        return new ParametricSampler<>(data, true);
    }

    private ParametricSampler(MolDataValues<T> data, boolean seed) {
        super(data, seed);
        gaussian = ZigguratSampler.NormalizedGaussian.of(rng);
    }

    /**
     * {@inheritDoc}
     *
     * <p>Returns a copy of the wrapped data with each observable independently
     * perturbed by Gaussian noise scaled by the per-field experimental error.
     */
    @Override
    public MolDataValues<T> sample() {
        MolDataValues<T> copy = data.copy();
        for (T dataValue : copy.getData()) {
            double[] originals = dataValue.getObservables();
            double[] errors = dataValue.getObservableErrors();
            double[] perturbed = new double[originals.length];
            for (int j = 0; j < perturbed.length; j++) {
                perturbed[j] = originals[j] + errors[j] * gaussian.sample();
            }
            dataValue.setObservables(perturbed);
        }
        copy.clearJValues();
        return copy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public MolDataValues<T> getOriginalData() {
        return data;
    }
}
