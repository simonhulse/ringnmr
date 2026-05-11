package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.List;

/**
 * Classifies the unique spectral-density frequencies of a deuterium dataset into
 * J(ω_D) and J(2ω_D) groups, replicating the field-list construction of
 * {@link DeuteriumMolDataValues#getNSpectralDensities()} and additionally tracking
 * the type of each non-zero entry.
 *
 * <p>Degeneracy rule: if ω_D of one field coincides (within 1 %) with 2ω_D of
 * another field, that frequency is placed in the J(ω_D) group.
 *
 * <p>The resulting {@link #wdIndices} and {@link #w2dIndices} arrays hold the
 * weight-vector indices (i.e., positions in the array returned by
 * {@link WeightSampler#sampleWeights()}) that belong to each group. Index 0 is
 * always J(0) and is not included in either group.
 */
class DeuteriumFrequencyGroups {

    final int[] wdIndices;
    final int[] w2dIndices;

    private DeuteriumFrequencyGroups(int[] wdIndices, int[] w2dIndices) {
        this.wdIndices  = wdIndices;
        this.w2dIndices = w2dIndices;
    }

    static DeuteriumFrequencyGroups of(List<DeuteriumDataValue> dataValues) {
        List<Double>  fieldList = new ArrayList<>();
        // isWd[i] is the type of fieldList[i+1]: true = wD, false = 2wD
        List<Boolean> isWd     = new ArrayList<>();
        fieldList.add(0.0);

        for (DeuteriumDataValue dv : dataValues) {
            double omega  = dv.getB0() * RelaxEquations.GAMMA_D / RelaxEquations.GAMMA_H * 2.0 * Math.PI;
            double omega2 = omega * 2.0;

            int omegaPos = findPos(fieldList, omega);
            if (omegaPos == -1) {
                fieldList.add(omega);
                isWd.add(true);
            } else {
                // Already present (was added as 2wD for a previous field) — degenerate → wD
                isWd.set(omegaPos - 1, true);
            }

            int omega2Pos = findPos(fieldList, omega2);
            if (omega2Pos == -1) {
                fieldList.add(omega2);
                isWd.add(false);
            }
            // If omega2 already exists, its classification was already resolved; leave it.
        }

        List<Integer> wdList  = new ArrayList<>();
        List<Integer> w2dList = new ArrayList<>();
        for (int i = 0; i < isWd.size(); i++) {
            (isWd.get(i) ? wdList : w2dList).add(i + 1);
        }

        return new DeuteriumFrequencyGroups(
            wdList.stream().mapToInt(Integer::intValue).toArray(),
            w2dList.stream().mapToInt(Integer::intValue).toArray()
        );
    }

    private static int findPos(List<Double> fieldList, double freq) {
        for (int i = 1; i < fieldList.size(); i++) {
            if (Math.abs((freq - fieldList.get(i)) / fieldList.get(i)) < 0.01) return i;
        }
        return -1;
    }
}
