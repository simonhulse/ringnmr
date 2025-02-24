package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;


public abstract class RandomSampler extends Sampler {
    int minSamples;
    int maxSamples;

    final public int getNSamples() {
        return ThreadLocalRandom.current().nextInt(minSamples, maxSamples + 1);
    }

    final public Double sampleUnwrapped() {
        if (minSamples != 1 || maxSamples != 1) {
            return null;
        } else {
            return sample().get(0);
        }
    }
}
