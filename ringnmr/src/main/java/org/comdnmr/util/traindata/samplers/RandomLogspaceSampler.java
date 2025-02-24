package org.comdnmr.util.traindata.samplers;

import java.util.*;


public class RandomLogspaceSampler extends RandomSampler {
    double minFirst;
    double maxFirst;
    double scale;

    public RandomLogspaceSampler(
        double minFirst,
        double maxFirst,
        double scale,
        int minSamples,
        int maxSamples
    ) {
        this.minFirst = minFirst;
        this.maxFirst = maxFirst;
        this.scale = scale;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    public RandomLogspaceSampler(
        String name,
        double minFirst,
        double maxFirst,
        double scale,
        int minSamples,
        int maxSamples
    ) {
        this(minFirst, maxFirst, scale, minSamples, maxSamples);
        this.name = name;
    }

    @Override
    public List<Double> sample() {
        int nSamples = getNSamples();
        var firstSampler = new RandomDoubleSampler(minFirst, maxFirst, 1, 1);
        double curr = firstSampler.sampleUnwrapped();
        List<Double> result = new ArrayList<Double>();

        for (int i = 0; i < nSamples; i++) {
            result.add(curr);
            curr *= scale;
        }

        return result;
    }
}
