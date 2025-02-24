package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class ExplicitSampler extends Sampler {
    List<Double> value;

    public ExplicitSampler(List<Double> value) {
        this.value = value;
    }

    public ExplicitSampler(String name, List<Double> value) {
        this(value);
        this.name = name;
    }

    @Override
    public List<Double> sample() {
        return value;
    }

    public Double sampleUnwrapped() {
        if (value.size() != 1) {
            return null;
        }
        return value.get(0);
    }
}
