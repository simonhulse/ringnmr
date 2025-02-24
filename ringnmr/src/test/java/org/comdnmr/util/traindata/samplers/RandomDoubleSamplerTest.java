package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.stream.Collectors;
import org.junit.Assert;
import org.junit.Test;


public class RandomDoubleSamplerTest {
    @Test
    public void testSampler() {
        // Generate 10 samplers and check their output is always appropriate
        // given the constructor arguments
        for (int i = 0; i < 10; i++) {
            Sampler sampler = new RandomDoubleSampler(0.0, 1.0, 2, 4);
            List<Double> sample = sampler.sample();
            Assert.assertTrue(sample.size() >= 2);
            Assert.assertTrue(sample.size() <= 4);
            Assert.assertEquals(sample, sample.stream().sorted().collect(Collectors.toList()));
            for (int n = 0; n < sample.size(); n++) {
                Assert.assertTrue(sample.get(n) >= 0.0);
                Assert.assertTrue(sample.get(n) <= 1.0);
            }
        }
    }
}
