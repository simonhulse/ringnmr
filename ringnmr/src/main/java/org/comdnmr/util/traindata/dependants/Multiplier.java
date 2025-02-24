package org.comdnmr.util.traindata.dependants;

import org.comdnmr.util.traindata.samplers.Sampler;

public class Multiplier extends DependantGenerator {
    Sampler sampler;
    double factor;

    public Multiplier(String name, Sampler sampler) {
        this.name = name;
        this.sampler = sampler;
        this.factor = 0.0;
    }

    public void updateState() {
        factor = sampler.sample().get(0);
    }

    public double fetch(double value) {
        return value * factor;
    }
}
