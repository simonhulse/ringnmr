package org.comdnmr.util.traindata.samplers;

import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

public abstract class Sampler {
    public String name;

    final public String getName() { return name; }
    abstract public List<Double> sample();
}
