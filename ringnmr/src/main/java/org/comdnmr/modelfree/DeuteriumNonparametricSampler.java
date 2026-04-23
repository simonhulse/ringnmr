package org.comdnmr.modelfree;

public class DeuteriumNonparametricSampler extends WeightSampler<DeuteriumDataValue> {

    public DeuteriumNonparametricSampler(MolDataValues<DeuteriumDataValue> data) {
        super(data);
    }

    @Override
    protected double[] sampleWeights() {
        int nJ = getNSpectralDensities();
        double[] weights = new double[nJ];
        weights[0] = 1.0;
        int nNonZero = nJ - 1;
        int i = 0;
        while (i < nNonZero) {
            int idx = 1 + rng.nextInt(nNonZero);
            if (weights[idx] < 2.0) {
                weights[idx] += 1.0;
                i++;
            }
        }
        return weights;
    }
}
