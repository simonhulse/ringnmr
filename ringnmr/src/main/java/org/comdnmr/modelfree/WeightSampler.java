package org.comdnmr.modelfree;

public abstract class WeightSampler extends BootstrapSampler {

    WeightSampler(MolDataValues data) { super(data); }

    protected abstract double[] sampleWeights();

    public MolDataValues sample() {
        double[] weights = sampleWeights();
        data.weight(weights);
        return data;
    }

    public MolDataValues getOriginalData() {
        int nJ = getNValues();
        double[] weights = new double[nJ];
        for (int j = 0; j < nJ; j++) weights[j] = 1.0;
        data.weight(weights);
        return data;
    }

}
