package org.comdnmr.util.traindata.samplers;

import java.util.*;


public class RandomSegmentedLinspaceSampler extends RandomSegmentedSampler {

    public RandomSegmentedLinspaceSampler(
        List<Double> minFirsts,
        List<Double> maxFirsts,
        List<Double> minLasts,
        List<Double> maxLasts,
        List<Integer> minSampless,
        List<Integer> maxSampless
    ) {
        samplers = new ArrayList<>();
        for (int i = 0; i < minFirsts.size(); i++) {
            double minFirst = minFirsts.get(i);
            double maxFirst = maxFirsts.get(i);
            double minLast = minLasts.get(i);
            double maxLast = maxLasts.get(i);
            int minSamples = minSampless.get(i);
            int maxSamples = maxSampless.get(i);

            var sampler = new RandomLinspaceSampler(minFirst, maxFirst, minLast, maxLast, minSamples, maxSamples);
            samplers.add(sampler);
        }
    }

    public RandomSegmentedLinspaceSampler(
        String name,
        List<Double> minFirsts,
        List<Double> maxFirsts,
        List<Double> minLasts,
        List<Double> maxLasts,
        List<Integer> minSampless,
        List<Integer> maxSampless
    ) {
        this(minFirsts, maxFirsts, minLasts, maxLasts, minSampless, maxSampless);
        this.name = name;
    }
}
