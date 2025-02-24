package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.stream.Collectors;
import org.junit.Assert;
import org.junit.Test;


public class RandomLinspaceSamplerTest {
    @Test
    public void testSampler() {
        // Generate 10 samplers and check their output is always appropriate
        // given the constructor arguments
        for (int i = 0; i < 10; i++) {
            Sampler sampler = new RandomLinspaceSampler(0.0, 1.0, 5.0, 10.0, 3, 6);
            List<Double> sample = sampler.sample();
            Assert.assertTrue(sample.size() >= 3);
            Assert.assertTrue(sample.size() <= 6);
            Assert.assertEquals(sample, sample.stream().sorted().collect(Collectors.toList()));
            double first = sample.get(0);
            double second = sample.get(1);
            double difference = second - first;
            for (int n = 1; n < sample.size() - 1; n++) {
                Assert.assertEquals(sample.get(n + 1) - sample.get(n), difference, 1.0e-10);
            }
        }
    }
}
