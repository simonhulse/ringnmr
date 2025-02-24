package org.comdnmr.util.traindata.samplers;

import java.util.*;


public class RandomChoiceSampler extends RandomSampler {
    List<Double> options;

    public RandomChoiceSampler(
        List<Double> options,
        int minSamples,
        int maxSamples
    ) {
        this.options = options;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    public RandomChoiceSampler(
        String name,
        List<Double> options,
        int minSamples,
        int maxSamples
    ) {
        this(options, minSamples, maxSamples);
        this.name = name;
    }

    @Override
    public List<Double> sample() {
        int nSamples = getNSamples();
        Collections.shuffle(options);
        List<Double> value = new ArrayList<>();

        for (int i = 0; i < nSamples; i++) {
            value.add(options.get(i));
        }

        Collections.sort(value);
        return value;
    }
}
