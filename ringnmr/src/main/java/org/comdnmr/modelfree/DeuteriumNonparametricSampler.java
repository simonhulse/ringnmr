package org.comdnmr.modelfree;

import java.util.*;

/**
 * Nonparametric bootstrap sampler for deuterium relaxation data.
 *
 * <p>For {@code nJ = getNSpectralDensities()} spectral densities, every valid weight
 * vector has {@code weights[0] = 1.0}, all other elements in {0.0, 1.0, 2.0}, and
 * total sum equal to {@code nJ}. The complete set of such vectors is enumerated at
 * construction time, shuffled, and iterated in order. Calling {@link #sample()} more
 * than the number of distinct weight vectors throws {@link java.util.NoSuchElementException}.
 *
 * <p>The total number of distinct weight vectors is:
 * \[ N(n_J) = \sum_{k=0}^{\lfloor (n_J-1)/2 \rfloor} \binom{n_J-1}{2k} \binom{2k}{k} \]
 * e.g. \(N(3)=3\), \(N(5)=19\), \(N(7)=141\), \(N(8)=393\), \(N(9)=1107\).
 */
public class DeuteriumNonparametricSampler extends WeightSampler<DeuteriumDataValue> {

    private final Iterator<int[]> iterator;
    private final int nBootstraps;

    public DeuteriumNonparametricSampler(MolDataValues<DeuteriumDataValue> data) {
        super(data);
        List<int[]> allWeights = enumerate(getNSpectralDensities());
        Collections.shuffle(allWeights, new Random(rng.nextLong()));
        nBootstraps = allWeights.size();
        iterator = allWeights.iterator();
    }

    public int getNBootstraps() { return nBootstraps; }

    @Override
    protected double[] sampleWeights() {
        if (!iterator.hasNext()) {
            throw new NoSuchElementException(String.format(
                "Maximum number of samples exceeded. " +
                "This sampler cannot be sampled more than %d times.",
                nBootstraps));
        }
        int[] src = iterator.next();
        double[] weights = new double[src.length + 1];
        weights[0] = 1.0;
        for (int i = 0; i < src.length; i++) weights[i + 1] = src[i];
        return weights;
    }

    private static List<int[]> enumerate(int nJ) {
        List<int[]> result = new ArrayList<>();
        fill(new int[nJ - 1], 0, nJ - 1, result);
        return result;
    }

    private static void fill(int[] suffix, int pos, int remaining, List<int[]> result) {
        if (pos == suffix.length) {
            if (remaining == 0) result.add(suffix.clone());
            return;
        }
        int slotsLeft = suffix.length - pos - 1;
        int lo = Math.max(0, remaining - 2 * slotsLeft);
        int hi = Math.min(2, remaining);
        for (int v = lo; v <= hi; v++) {
            suffix[pos] = v;
            fill(suffix, pos + 1, remaining - v, result);
        }
    }
}
