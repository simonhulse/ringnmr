package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.stream.Collectors;
import org.junit.Assert;
import org.junit.Test;


public class RandomChoiceSamplerTest {
    @Test
    public void testSampler() {
        // Generate 10 samplers and check their output is always appropriate
        // given the constructor arguments
        List<Double> options = Arrays.asList(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
        for (int i = 0; i < 10; i++) {
            Sampler sampler = new RandomChoiceSampler(options, 3, 5);
            List<Double> sample = sampler.sample();
            Assert.assertTrue(sample.size() >= 3);
            Assert.assertTrue(sample.size() <= 5);
            Assert.assertEquals(sample, sample.stream().sorted().collect(Collectors.toList()));
            for (int n = 0; n < sample.size(); n++) {
                Assert.assertTrue(options.contains(sample.get(n)));
            }
        }
    }
}
