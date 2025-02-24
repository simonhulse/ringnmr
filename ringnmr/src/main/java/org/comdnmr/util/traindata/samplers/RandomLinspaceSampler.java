package org.comdnmr.util.traindata.samplers;

import java.util.*;


public class RandomLinspaceSampler extends RandomSampler {
    double minFirst;
    double maxFirst;
    double minLast;
    double maxLast;

    public RandomLinspaceSampler(
        double minFirst,
        double maxFirst,
        double minLast,
        double maxLast,
        int minSamples,
        int maxSamples
    ) {
        this.minFirst = minFirst;
        this.maxFirst = maxFirst;
        this.minLast = minLast;
        this.maxLast = maxLast;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    public RandomLinspaceSampler(
        String name,
        double minFirst,
        double maxFirst,
        double minLast,
        double maxLast,
        int minSamples,
        int maxSamples
    ) {
        this(minFirst, maxFirst, minLast, maxLast, minSamples, maxSamples);
        this.name = name;
    }

    @Override
    public List<Double> sample() {
        RandomDoubleSampler firstSampler = new RandomDoubleSampler(minFirst, maxFirst, 1, 1);
        double first = firstSampler.sampleUnwrapped();

        RandomDoubleSampler lastSampler = new RandomDoubleSampler(minLast, maxLast, 1, 1);
        double last = lastSampler.sampleUnwrapped();

        int nSamples = getNSamples();

        double incr = (last - first) / (nSamples - 1);
        double curr = first;

        List<Double> value = new ArrayList<Double>();
        for (int i = 0; i < nSamples; i++) {
            value.add(curr);
            curr += incr;
        }

        return value;
    }
}
