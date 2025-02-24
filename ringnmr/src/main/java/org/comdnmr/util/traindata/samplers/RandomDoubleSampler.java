package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;


public class RandomDoubleSampler extends RandomSampler {
    double min;
    double max;

    public RandomDoubleSampler(double min, double max, int minSamples, int maxSamples) {
        this.min = min;
        this.max = max;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    public RandomDoubleSampler(String name, double min, double max, int minSamples, int maxSamples) {
        this(min, max, minSamples, maxSamples);
        this.name = name;
    }

    public List<Double> sample() {
        int nSamples = getNSamples();
        List<Double> value = new ArrayList<Double>();
        for (int i = 0; i < nSamples; i++) {
            value.add(ThreadLocalRandom.current().nextDouble(min, max));
        }
        Collections.sort(value);
        return value;
    }
}
