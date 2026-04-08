package org.comdnmr.modelfree;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.ObjectSampler;
import org.apache.commons.rng.simple.RandomSource;


public abstract class BootstrapSampler implements ObjectSampler<MolDataValues> {

    private static final int[] SEED = new int[] {196, 9, 0, 226};

    protected final MolDataValues data;
    protected final UniformRandomProvider rng;

    public BootstrapSampler(MolDataValues data) { this(data, false); }

    public BootstrapSampler(MolDataValues data, boolean seed) {
        this.data = data;
        this.rng = (seed) ?
            RandomSource.XO_SHI_RO_128_PP.create(SEED) :
            RandomSource.XO_SHI_RO_128_PP.create();
    }

    public int getNValues() {
        return data.getNValues();
    }

    public int getNFields() {
        return data.dataValues.size();
    }

    public abstract MolDataValues sample();

    public abstract MolDataValues getOriginalData();
}
