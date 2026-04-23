package org.comdnmr.modelfree;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.ObjectSampler;
import org.apache.commons.rng.simple.RandomSource;


/**
 * Abstract base class for bootstrap resampling strategies used in model-free analysis.
 *
 * <p>A {@code BootstrapSampler} wraps a {@link MolDataValues} object and produces
 * resampled copies of it for uncertainty quantification. Each call to {@link #sample()}
 * returns the <em>same</em> underlying {@code MolDataValues} reference with its state
 * mutated to reflect the current bootstrap replicate; callers must therefore process or
 * copy the result before drawing the next sample.
 *
 * <p>Two concrete sub-hierarchies exist:
 * <ul>
 *   <li>{@link WeightSampler} (and its subclasses {@link NonparametricSampler} and
 *       {@link BayesianSampler}) — resamples by assigning per-observation weights
 *       without modifying the measured values.</li>
 *   <li>{@link ParametricSampler} — resamples by drawing new measured values from
 *       Gaussian distributions centred on the originals, using the reported
 *       experimental errors as standard deviations.</li>
 * </ul>
 *
 * <p>The RNG is an XoShiRo128++ generator from Apache Commons RNG. Passing
 * {@code seed = true} to the two-argument constructor fixes the seed to a
 * hard-coded value for reproducible results (e.g., in tests).
 */
public abstract class BootstrapSampler<T extends RelaxDataValue> implements ObjectSampler<MolDataValues<T>> {

    private static final int[] SEED = new int[] {196, 9, 0, 226};

    protected final MolDataValues<T> data;
    protected final UniformRandomProvider rng;

    /**
     * Constructs a {@code BootstrapSampler} with a randomly seeded RNG.
     *
     * @param data the relaxation data to be resampled
     */
    public BootstrapSampler(MolDataValues<T> data) { this(data, false); }

    /**
     * Constructs a {@code BootstrapSampler}, optionally fixing the RNG seed.
     *
     * @param data the relaxation data to be resampled
     * @param seed if {@code true}, initialises the RNG with a fixed seed for
     *             reproducibility; if {@code false}, uses a random seed
     */
    public BootstrapSampler(MolDataValues<T> data, boolean seed) {
        this.data = data;
        this.rng = (seed) ?
            RandomSource.XO_SHI_RO_128_PP.create(SEED) :
            RandomSource.XO_SHI_RO_128_PP.create();
    }

    /**
     * Returns the total number of weighted spectral-density observations in the wrapped
     * data (i.e. {@code 3 × nFields} for R1/R2/NOE data).
     *
     * @return the number of J(ω) values
     */
    public int getNSpectralDensities() {
        return data.getNSpectralDensities();
    }

    /**
     * Returns the number of static magnetic fields at which data were acquired, i.e.
     * the number of {@link RelaxDataValue} entries in the wrapped {@link MolDataValues}.
     *
     * @return the number of relaxation data entries
     */
    public int getNFields() {
        return data.getData().size();
    }

    /**
     * Draws the next bootstrap replicate by mutating the wrapped {@link MolDataValues}
     * in-place and returning it.
     *
     * <p><strong>Important:</strong> the returned reference is the same object on every
     * call. Callers that need to retain the state of a replicate must copy the result
     * before drawing the next sample.
     *
     * @return the wrapped {@code MolDataValues} in its resampled state
     */
    public abstract MolDataValues<T> sample();

    /**
     * Restores the wrapped {@link MolDataValues} to its original, unperturbed state and
     * returns it.
     *
     * @return the wrapped {@code MolDataValues} with the original data restored
     */
    public abstract MolDataValues<T> getOriginalData();
}
