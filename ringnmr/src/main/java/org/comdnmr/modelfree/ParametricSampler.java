package org.comdnmr.modelfree;

import java.util.List;
import java.util.stream.Collectors;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratSampler;

/**
 * Parametric bootstrap sampler that generates synthetic observations by adding
 * Gaussian noise to the original measured values.
 *
 * <p>On each call to {@link #sample()}, every relaxation observable (R1, R2,
 * and either NOE for {@link R1R2NOEDataValue} or rQ/rAP for
 * {@link DeuteriumDataValue}) is independently perturbed by a zero-mean Gaussian
 * with standard deviation equal to the corresponding reported experimental error.
 * The original values are preserved internally and restored by
 * {@link #getOriginalData()}.
 *
 * <p>Normal deviates are drawn using the Ziggurat algorithm via Apache Commons RNG.
 *
 * <p><strong>Note:</strong> the underlying {@link MolDataValues} object is mutated
 * in-place. After calling {@link #sample()}, it holds perturbed values until the
 * next call to {@link #sample()} or {@link #getOriginalData()}.
 */
public class ParametricSampler<T extends RelaxDataValue> extends BootstrapSampler<T> {

    private final NormalizedGaussianSampler gaussian;
    private final List<double[]> originalObservables;

    /**
     * Constructs a {@code ParametricSampler} and snapshots the original observed
     * values from {@code data} for later restoration.
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
        originalObservables = data.getData().stream()
            .map(T::getObservables)
            .collect(Collectors.toList());
    }

    /**
     * {@inheritDoc}
     *
     * <p>Perturbs each observable by independent Gaussian noise scaled by the
     * per-field experimental error. Clears the cached J-values so they are
     * recomputed on next access.
     */
    public MolDataValues<T> sample() {
        List<T> dataValues = data.getData();
        for (int i = 0; i < getNFields(); i++) {
            T dataValue = dataValues.get(i);
            double[] originals = originalObservables.get(i);
            double[] errors = dataValue.getObservableErrors();
            double[] perturbed = new double[originals.length];
            for (int j = 0; j < perturbed.length; j++) {
                perturbed[j] = originals[j] + errors[j] * gaussian.sample();
            }
            dataValue.setObservables(perturbed);
        }
        data.clearJValues();
        return data;
    }

    /**
     * {@inheritDoc}
     *
     * <p>Restores every observable to the value snapshotted at construction and
     * clears the cached J-values.
     */
    public MolDataValues<T> getOriginalData() {
        List<T> dataValues = data.getData();
        for (int i = 0; i < getNFields(); i++) {
            dataValues.get(i).setObservables(originalObservables.get(i));
        }
        data.clearJValues();
        return data;
    }
}
