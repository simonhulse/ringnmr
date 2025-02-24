package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.stream.Collectors;
import org.junit.Assert;
import org.junit.Test;


public class RandomLogspaceSamplerTest {
    @Test
    public void testSampler() {
        // Generate 10 samplers and check their output is always appropriate
        // given the constructor arguments
        double scale = 2.0;
        for (int i = 0; i < 10; i++) {
            Sampler sampler = new RandomLogspaceSampler(0.0, 1.0, scale, 3, 6);
            List<Double> sample = sampler.sample();
            Assert.assertTrue(sample.size() >= 3);
            Assert.assertTrue(sample.size() <= 6);
            Assert.assertEquals(sample, sample.stream().sorted().collect(Collectors.toList()));
            System.out.println(String.format("sample: %s", sample));
            for (int n = 0; n < sample.size() - 1; n++) {
                Assert.assertEquals(sample.get(n + 1) / sample.get(n), scale, 1.0e-10);
            }
        }
    }
}
