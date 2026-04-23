package org.comdnmr.modelfree;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Nonparametric bootstrap sampler based on the spectral-density selection patterns
 * described by Crawley and Palmer (2021).
 *
 * <p>For {@code nFields} static magnetic fields, a discrete table ({@link #SELECTION_MAP})
 * enumerates every structurally valid resampling of the J(ω) values. A single bootstrap
 * replicate is formed by choosing one row from this table independently for each of the
 * three canonical spectral-density frequency classes — J(0), J(ω_N), J(0.87ω_H) —
 * giving {@code nSel³} distinct weight vectors in total.
 *
 * <p>The full index space is shuffled once at construction time and then iterated in
 * random order. Calling {@link #sample()} more than {@code nSel³} times throws
 * {@link java.util.NoSuchElementException}.
 *
 * <p>Each weight vector element encodes how many times the corresponding J(ω)
 * observation is selected across the three frequency draws; selected twice → weight 2.0,
 * not selected → weight 0.0.
 */
public class NonparametricSampler<T extends RelaxDataValue> extends WeightSampler<T> {

    // Iterator of randomly ordered ints which specifies the ordering of possible weight vectors.
    private final Iterator<Integer> iterator;

    // Static map containing permitted combinations of weight assignments for a
    // given frequency class. See Table 1 in Crawley and Palmer's paper.
    private static final Map<Integer, int[][]> SELECTION_MAP;

    static {
        SELECTION_MAP = new HashMap<>();
        SELECTION_MAP.put(
            2,
            new int[][]{
                {0, 0},
                {0, 1},
                {1, 1}
            }
        );
        SELECTION_MAP.put(
            3,
            new int[][]{
                {0, 1, 2},
                {1, 1, 2},
                {0, 0, 2},
                {2, 1, 2},
                {0, 1, 0},
                {0, 2, 2},
                {0, 1, 1}
            }
        );
        SELECTION_MAP.put(
            4,
            new int[][]{
                {0, 1, 2, 3},
                {0, 0, 2, 3},
                {0, 1, 0, 3},
                {0, 1, 2, 0},
                {1, 1, 2, 3},
                {0, 1, 1, 3},
                {0, 1, 2, 1},
                {2, 1, 2, 3},
                {0, 2, 2, 3},
                {0, 1, 2, 2},
                {3, 1, 2, 3},
                {0, 3, 2, 3},
                {0, 1, 3, 3},
                {0, 0, 1, 1},
                {0, 0, 2, 2},
                {0, 0, 3, 3},
                {1, 1, 2, 2},
                {1, 1, 3, 3},
                {2, 2, 3, 3}
            }
        );
    }

    /**
     * Constructs a {@code NonparametricSampler} for the given relaxation data and
     * pre-shuffles the full bootstrap index space.
     *
     * <p>The number of static fields is inferred from {@code data.dataValues.size()};
     * it must be 2, 3, or 4 (matching the keys in {@link #SELECTION_MAP}).
     *
     * @param data the relaxation data to be resampled
     */
    public NonparametricSampler(MolDataValues<T> data) {
        super(data);
        iterator = generateIterator();
    }

    /**
     * Generates a randomized iterator over bootstrap sample indices.
     *
     * @return an iterator over the randomized bootstrap indices
     */
    private Iterator<Integer> generateIterator() {
        // Create a list of integers from 0 to N-1, shuffle it, and create an iterator
        List<Integer> randomNumbers = IntStream.range(0, getNBootstraps())
            .boxed()
            .collect(Collectors.toList());
        Collections.shuffle(randomNumbers, new Random(rng.nextLong()));
        return randomNumbers.iterator();
    }

    /**
     * Returns the number of selections based on the number of static fields.
     *
     * @return the number of selections
     */
    public int getNSelections() {
        return getSelections().length;
    }

    /**
     * Returns the selection patterns, given `nExp()` static fields.
     *
     * @return a 2D array of selection patterns
     */
    private int[][] getSelections() {
        return SELECTION_MAP.get(getNFields());
    }

    /**
     * Retrieves a specific row of selections based on the provided index.
     *
     * @param index the index of the selection row
     * @return an array representing the selection row
     */
    private int[] getSelectionsRow(int index) {
        return getSelections()[index];
    }

    /**
     * Calculates the total number of bootstrap samples.
     *
     * @return the total number of bootstraps
     */
    public int getNBootstraps() {
        int nSel = getNSelections();
        return nSel * nSel * nSel;
    }

    /**
     * Returns the next weight vector from the pre-shuffled bootstrap index space.
     *
     * <p>Each element of the returned array is the number of times the corresponding
     * J(ω) observation was selected (0, 1, 2, or 3) across the three frequency draws.
     *
     * @return per-observation bootstrap weights
     * @throws NoSuchElementException if the sampler has been exhausted
     */
    protected double[] sampleWeights() {
        if (!iterator.hasNext()) {
            throw new NoSuchElementException(
                String.format(
                    "Maximum number of samples exceeded. " +
                    "This sampler cannot be sampled more than %d (%d^3) times.",
                    getNBootstraps(),
                    getNSelections()
                )
            );
        }
        int index = iterator.next();

        int nSel = getNSelections();
        int nSelSq = nSel * nSel;

        // J(0) pointer
        int p0 = index / nSelSq;
        int r0 = index % nSelSq;
        // J(wN) pointer
        int pN = r0 / nSel;
        // J(0.87wH) pointer
        int pH = r0 % nSel;

        int[][] selections = {
            getSelectionsRow(p0),
            getSelectionsRow(pN),
            getSelectionsRow(pH)
        };

        double[] weights = new double[getNSpectralDensities()];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < getNFields(); j++) {
                weights[selections[i][j] * 3 + i] += 1.0;
            }
        }
        return weights;
    }
}
