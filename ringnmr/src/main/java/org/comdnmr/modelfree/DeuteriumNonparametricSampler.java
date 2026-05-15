package org.comdnmr.modelfree;

import java.util.*;

/**
 * Nonparametric bootstrap sampler for deuterium relaxation data.
 *
 * <p>The J(ω) observations are divided into three groups: J(0) receives a fixed
 * weight of 1, while J(ω_D) and J(2ω_D) observations are each resampled
 * independently. For a group of size M, every valid weight vector — integer
 * elements in {0, 1, 2} summing to M — is enumerated. A bootstrap replicate is
 * formed by choosing one vector from each group; the full sample space is therefore
 * the Cartesian product of the two group enumerations.
 *
 * <p>Degenerate frequencies (where ω_D of one field coincides with 2ω_D of another)
 * are placed in the J(ω_D) group; see {@link DeuteriumFrequencyGroups}.
 *
 * <p>The Cartesian product is shuffled once at construction time and then iterated
 * in random order. Calling {@link #sample()} more than {@link #getNBootstraps()} times
 * throws {@link java.util.NoSuchElementException}.
 *
 * <p>Bootstrap counts for representative field configurations:
 * <ul>
 *   <li>1 field: 1</li>
 *   <li>2 fields, no degeneracy: 9</li>
 *   <li>2 fields, B = 2A: 3</li>
 *   <li>3 fields, no degeneracy: 49</li>
 *   <li>3 fields, one pair: 21</li>
 *   <li>4 fields, no degeneracy: 361</li>
 *   <li>4 fields, one pair: 133</li>
 *   <li>4 fields, two pairs / chain-of-3: 57</li>
 * </ul>
 */
public class DeuteriumNonparametricSampler extends WeightSampler<DeuteriumDataValue> {

    private final DeuteriumFrequencyGroups groups;
    private final List<int[]> wdVectors;
    private final List<int[]> w2dVectors;
    private final Iterator<int[]> iterator;
    private final int nBootstraps;

    public DeuteriumNonparametricSampler(MolDataValues<DeuteriumDataValue> data) { this(data, false); }

    // Added for use in the regularization paper; not used within RING.
    public static DeuteriumNonparametricSampler withFixedSeed(MolDataValues<DeuteriumDataValue> data) {
        return new DeuteriumNonparametricSampler(data, true);
    }

    private DeuteriumNonparametricSampler(MolDataValues<DeuteriumDataValue> data, boolean seed) {
        super(data, seed);
        groups     = DeuteriumFrequencyGroups.of(data.getData());
        wdVectors  = enumerate(groups.wdIndices.length);
        w2dVectors = enumerate(groups.w2dIndices.length);

        List<int[]> pairs = new ArrayList<>(wdVectors.size() * w2dVectors.size());
        for (int i = 0; i < wdVectors.size(); i++)
            for (int j = 0; j < w2dVectors.size(); j++)
                pairs.add(new int[]{i, j});
        Collections.shuffle(pairs, new Random(rng.nextLong()));

        nBootstraps = pairs.size();
        iterator    = pairs.iterator();
    }

    public int getNBootstraps() { return nBootstraps; }

    @Override
    public MolDataValues<DeuteriumDataValue> sample() {
        if (!iterator.hasNext()) {
            throw new NoSuchElementException(String.format(
                "Maximum number of samples exceeded. " +
                "This sampler cannot be sampled more than %d times.",
                nBootstraps));
        }
        int[] pair       = iterator.next();
        double[] weights = new double[getNSpectralDensities()];
        weights[0] = 1.0;
        int[] wd  = wdVectors.get(pair[0]);
        int[] w2d = w2dVectors.get(pair[1]);
        for (int i = 0; i < groups.wdIndices.length;  i++) weights[groups.wdIndices[i]]  = wd[i];
        for (int i = 0; i < groups.w2dIndices.length; i++) weights[groups.w2dIndices[i]] = w2d[i];
        MolDataValues<DeuteriumDataValue> copy = data.copy();
        copy.setWeights(weights);
        return copy;
    }

    /**
     * Enumerates all integer arrays of length {@code size} with elements in {0, 1, 2}
     * whose elements sum to {@code size}.
     */
    private static List<int[]> enumerate(int size) {
        List<int[]> result = new ArrayList<>();
        if (size == 0) {
            result.add(new int[0]);
            return result;
        }
        fill(new int[size], 0, size, result);
        return result;
    }

    private static void fill(int[] vec, int pos, int remaining, List<int[]> result) {
        if (pos == vec.length) {
            if (remaining == 0) result.add(vec.clone());
            return;
        }
        int slotsLeft = vec.length - pos - 1;
        int lo = Math.max(0, remaining - 2 * slotsLeft);
        int hi = Math.min(2, remaining);
        for (int v = lo; v <= hi; v++) {
            vec[pos] = v;
            fill(vec, pos + 1, remaining - v, result);
        }
    }
}
