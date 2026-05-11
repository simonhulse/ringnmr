package org.comdnmr.modelfree;

import java.util.*;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.stat.inference.TTest;

import org.comdnmr.data.DynamicsSource;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;


public class BootstrapSamplerTest {

    private static final double DELTA = 1.0e-12;
    private static final double ALPHA = 0.01;

    private static final Map<String, List<Double>> ASP_3_DATA = new HashMap<>();
    static {
        ASP_3_DATA.put("B0", List.of(600.133, 700.133, 800.234));
        ASP_3_DATA.put("R1", List.of(1.111043, 1.06656043, 1.091134));
        ASP_3_DATA.put("R1err", List.of(0.01628869, 0.0330645, 0.07346691));
        ASP_3_DATA.put("R2", List.of(3.05152, 1.47590564, 3.587807));
        ASP_3_DATA.put("R2err", List.of(0.38221198, 0.10852504, 0.04461263));
        ASP_3_DATA.put("NOE", List.of(-0.758704, -0.4115185, -0.175125));
        ASP_3_DATA.put("NOEerr", List.of(0.00120198, 0.00462507, 0.00075348));
    };

    private R1R2NOEMolDataValues data;

    @Before
    public void setUp() {
        int nFields = ASP_3_DATA.get("R1").size();
        data = new R1R2NOEMolDataValues("1.H", new double[3], new DynamicsSource(true, true, true, true));
        for (int i = 0; i < nFields; i++) {
            double r1 = ASP_3_DATA.get("R1").get(i);
            double r1Err = ASP_3_DATA.get("R1err").get(i);
            double r2 = ASP_3_DATA.get("R2").get(i);
            double r2Err = ASP_3_DATA.get("R2err").get(i);
            double noe = ASP_3_DATA.get("NOE").get(i);
            double noeErr = ASP_3_DATA.get("NOEerr").get(i);
            double b0 = ASP_3_DATA.get("B0").get(i);
            RelaxEquations relaxEquations = new RelaxEquations(b0, "H", "N");
            R1R2NOEDataValue dataValue = new R1R2NOEDataValue(data, r1, r1Err, r2, r2Err, noe, noeErr, relaxEquations);
            data.addData(dataValue);
        }
    }

    private double computeMean(double[] xs) {
        double mean = 0.0;
        for (double x : xs) {
            mean += x;
        }
        mean /= xs.length;
        return mean;
    }

    private double computeStdev(double[] xs) {
        double mean = computeMean(xs);
        double variance = 0.0;
        for (double x : xs) {
            variance += Math.pow(x - mean, 2.0);
        }
        variance /= xs.length - 1;
        double stdev = Math.sqrt(variance);
        return stdev;
    }

    private void testSample(String description, double[] sample, double expectedMean, double expectedStdev) {
        TTest tTest = new TTest();
        double pValueMean = tTest.tTest(expectedMean, sample);
        assertTrue(
            String.format(
                "%s:%nSample mean is significantly different for the expected mean.%np-value: %.3f",
                description,
                pValueMean
            ),
            pValueMean > ALPHA
        );

        // Test standard deviation with chi-squared
        int n = sample.length;
        double sampleStdev = computeStdev(sample);
        double chiSquared = ((n - 1) * Math.pow(sampleStdev, 2)) / Math.pow(expectedStdev, 2);

        ChiSquaredDistribution distribution = new ChiSquaredDistribution(n - 1);
        double pValueStdev = 1.0 - distribution.cumulativeProbability(chiSquared);
        assertTrue(
            String.format(
                "%s:%nSample stdev of R1 samples is significantly different for the expected stdev.%np-value: %.3f",
                description,
                pValueStdev
            ),
            pValueStdev > ALPHA
        );
    }

    @Test
    public void testParametricSampler() {
        ParametricSampler<R1R2NOEDataValue> sampler = new ParametricSampler<>(data);

        // Check that values of R1 sampled resemble normal distribution with
        // expected mean and stdev
        // Using t-test for mean and chi-squared test for stdev
        int nSamples = 50;
        double[] r1Samples = new double[nSamples];
        for (int i = 0; i < nSamples; i++) {
            MolDataValues<R1R2NOEDataValue> bootstrapData = sampler.sample();
            r1Samples[i] = bootstrapData.getData().get(0).R1;
        }

        double expectedMean = ASP_3_DATA.get("R1").get(0);
        double expectedStdev = ASP_3_DATA.get("R1err").get(0);

        testSample(
            String.format("ASP-3 | R1 | B0 = %.2f MHz", ASP_3_DATA.get("B0").get(0)),
            r1Samples,
            expectedMean,
            expectedStdev
        );

        // Check that original data can be recovered
        MolDataValues<R1R2NOEDataValue> originalData = sampler.getOriginalData();
        int nFields = originalData.getData().size();
        for (int i = 0; i < nFields; i++) {
            R1R2NOEDataValue value = originalData.getData().get(i);
            assertEquals(ASP_3_DATA.get("R1").get(i), value.R1, DELTA);
            assertEquals(ASP_3_DATA.get("R1err").get(i), value.R1err, DELTA);
            assertEquals(ASP_3_DATA.get("R2").get(i), value.R2, DELTA);
            assertEquals(ASP_3_DATA.get("R2err").get(i), value.R2err, DELTA);
            assertEquals(ASP_3_DATA.get("NOE").get(i), value.NOE, DELTA);
            assertEquals(ASP_3_DATA.get("NOEerr").get(i), value.NOEerr, DELTA);
        }
    }

    private boolean checkNonparametricWeightVector(double[] weights) {
        for (int j = 0; j < weights.length; j++) {
            double weight = weights[j];
            if (
                !(
                    Math.abs(weight - 0.0) < DELTA ||
                    Math.abs(weight - 1.0) < DELTA ||
                    Math.abs(weight - 2.0) < DELTA
                )
            ) {
                return false;
            }
        }
        return checkWeightVectorHasCorrectSum(weights);
    }

    private boolean checkWeightVectorHasCorrectSum(double weights[]) {
        double sum = 0.0;
        for (int j = 0; j < weights.length; j++) {
            sum += weights[j];
        }
        return Math.abs(sum - weights.length) < DELTA;
    }

    private boolean checkVectorOfOnes(double[] weights) {
        for (int j = 0; j < weights.length; j++) {
            if (!(Math.abs(weights[j] - 1.0) < DELTA)) return false;
        }
        return true;
    }

    private int generateWeightHash(double[] weights) {
        // `weights` array is treated like a ternary number, with hash computed
        // as decimal representation
        int hash = 0;
        int nWeights = weights.length;
        for (int i = 0; i < nWeights; i++) {
            int n = nWeights - i - 1;
            hash += (int) weights[i] * Math.pow(3, n);
        }
        return hash;
    }

    @Test
    public void testAmideNonparametricSampler() {
        AmideNonparametricSampler sampler = new AmideNonparametricSampler(data);
        int nIterations = 343;
        Set<Integer> weightSet = new HashSet<>();
        for (int i = 0; i < nIterations; i++) {
            MolDataValues<R1R2NOEDataValue> sample = sampler.sample();
            double[] weights = sample.getWeights();

            // Check that all weights are 0.0, 1.0 or 2.0
            assertTrue(
                String.format(
                    "Weight vector doesn't comprise solely 0.0, 1.0 and 2.0, and/or doesn't have a sum of %d",
                    weights.length
                ),
                checkNonparametricWeightVector(weights)
            );

            // Check the weight vector has not already been sampled
            int weightHash = generateWeightHash(weights);
            assertFalse(
                String.format("Weight vector %s already sampled", Arrays.toString(weights)),
                weightSet.contains(weightHash)
            );
            weightSet.add(weightHash);
        }

        // After 343 iterations, should have exhausted the sampler
        try {
            sampler.sample();
            fail(
                String.format(
                    "Expected nonparameteric sampler to be exhausted after %d iterations",
                    nIterations
                )
            );
        } catch (NoSuchElementException e) {
            assertEquals(
                "Maximum number of samples exceeded. This sampler cannot be sampled more than 343 (7^3) times.",
                e.getMessage()
            );
        }

        // Check original data has weight vector of ones
        MolDataValues<R1R2NOEDataValue> originalData = sampler.getOriginalData();
        assertTrue(
            "Original data should have a weight vector of ones",
            checkVectorOfOnes(originalData.getWeights())
        );
    }

    @Test
    public void testBayesianSampler() {
        AmideBayesianSampler sampler = new AmideBayesianSampler(data);

        int nSamples = 100;
        for (int i = 0; i < nSamples; i++) {
            MolDataValues<R1R2NOEDataValue> sample = sampler.sample();
            double[] weights = sample.getWeights();
            assertTrue(
                String.format("Expected sum of weights to equal %d", weights.length),
                checkWeightVectorHasCorrectSum(weights)
            );
        }

        // Check original data has weight vector of ones
        MolDataValues<R1R2NOEDataValue> originalData = sampler.getOriginalData();
        assertTrue(
            "Original data should have a weight vector of ones",
            checkVectorOfOnes(originalData.getWeights())
        );
    }

}
