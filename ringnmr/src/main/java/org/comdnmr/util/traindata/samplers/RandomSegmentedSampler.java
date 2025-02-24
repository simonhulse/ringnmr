package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;


abstract class RandomSegmentedSampler extends Sampler {
    List<RandomSampler> samplers;

    @Override
    final public List<Double> sample() {
        List<Double> value = new ArrayList<>();
        for (RandomSampler sampler : samplers) {
            List<Double> subValue = sampler.sample();
            value.addAll(subValue);
        }
        Collections.sort(value);
        return value;
    }
}
