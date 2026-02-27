/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.comdnmr.fit.calc;

import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.eqnfit.FitResult;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import org.comdnmr.eqnfit.CESTFitter;
import org.comdnmr.util.CoMDOptions;
import org.junit.Assert;
import org.junit.Test;
import org.nmrfx.datasets.Nuclei;

public class CESTFitTest {

    private double[] refValues;

    private List<Double> offset = Arrays.asList(-6220.353454107791, -6031.857894892402, -5843.362335677015, -5654.8667764616275, -5466.37121724624,
            -5277.875658030853, -5089.380098815464, -4900.884539600077, -4712.38898038469, -4523.893421169302,
            -4335.397861953915, -4146.9023027385265, -3958.4067435231395, -3769.9111843077517, -3581.4156250923643,
            -3392.9200658769764, -3204.424506661589, -3015.928947446201, -2827.4333882308138, -2638.9378290154264,
            -2450.4422698000385, -2261.946710584651, -2073.4511513692632, -1884.9555921538758, -1696.4600329384882,
            -1507.9644737231006, -1319.4689145077132, -1130.9733552923256, -942.4777960769379, -753.9822368615503,
            -565.4866776461628, -376.99111843077515, -188.49555921538757, 0.0, 188.49555921538757, 376.99111843077515,
            565.4866776461628, 753.9822368615503, 942.4777960769379, 1130.9733552923256, 1319.4689145077132,
            1507.9644737231006, 1696.4600329384882, 1884.9555921538758, 2073.4511513692632, 2261.946710584651,
            2450.4422698000385, 2638.9378290154264, 2827.4333882308138, 3015.928947446201, 3204.424506661589,
            3392.9200658769764, 3581.4156250923643, 3769.9111843077517, 3958.4067435231395, 4146.9023027385265,
            4335.397861953915, 4523.893421169302, 4712.38898038469, 4900.884539600077, 5089.380098815464,
            5277.875658030853, 5466.37121724624, 5654.8667764616275, 5843.362335677015, 6031.857894892402,
            6220.353454107791, -6283.185307179586, -6031.857894892402, -5780.53048260522, -5529.203070318036,
            -5277.875658030853, -5026.548245743669, -4775.220833456486, -4523.893421169302, -4272.566008882119,
            -4021.238596594935, -3769.9111843077517, -3518.583772020568, -3267.2563597333847, -3015.928947446201,
            -2764.601535159018, -2513.2741228718346, -2261.946710584651, -2010.6192982974676, -1759.291886010284,
            -1507.9644737231006, -1256.6370614359173, -1005.3096491487338, -753.9822368615503, -502.6548245743669,
            -251.32741228718345, 0.0, 251.32741228718345, 502.6548245743669, 753.9822368615503, 1005.3096491487338,
            1256.6370614359173, 1507.9644737231006, 1759.291886010284, 2010.6192982974676, 2261.946710584651,
            2513.2741228718346, 2764.601535159018, 3015.928947446201, 3267.2563597333847, 3518.583772020568,
            3769.9111843077517, 4021.238596594935, 4272.566008882119, 4523.893421169302, 4775.220833456486,
            5026.548245743669, 5277.875658030853, 5529.203070318036, 5780.53048260522, 6031.857894892402,
            6283.185307179586, -6283.185307179586, -5969.026041820607, -5654.8667764616275, -5340.707511102648,
            -5026.548245743669, -4712.38898038469, -4398.22971502571, -4084.070449666731, -3769.9111843077517,
            -3455.7519189487725, -3141.592653589793, -2827.4333882308138, -2513.2741228718346, -2199.114857512855,
            -1884.9555921538758, -1570.7963267948965, -1256.6370614359173, -942.4777960769379, -628.3185307179587,
            -314.1592653589793, 0.0, 314.1592653589793, 628.3185307179587, 942.4777960769379, 1256.6370614359173,
            1570.7963267948965, 1884.9555921538758, 2199.114857512855, 2513.2741228718346, 2827.4333882308138,
            3141.592653589793, 3455.7519189487725, 3769.9111843077517, 4084.070449666731, 4398.22971502571,
            4712.38898038469, 5026.548245743669, 5340.707511102648, 5654.8667764616275, 5969.026041820607,
            6283.185307179586);

    private List<Double> omega = Arrays.asList(111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 111.08671623, 111.08671623, 111.08671623,
            111.08671623, 111.08671623, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 175.30087007, 175.30087007,
            175.30087007, 175.30087007, 175.30087007, 302.91236366, 302.91236366,
            302.91236366, 302.91236366, 302.91236366, 302.91236366, 302.91236366,
            302.91236366, 302.91236366, 302.91236366, 302.91236366, 302.91236366,
            302.91236366, 302.91236366, 302.91236366, 302.91236366, 302.91236366,
            302.91236366, 302.91236366, 302.91236366, 302.91236366, 302.91236366,
            302.91236366, 302.91236366, 302.91236366, 302.91236366, 302.91236366,
            302.91236366, 302.91236366, 302.91236366, 302.91236366, 302.91236366,
            302.91236366, 302.91236366, 302.91236366, 302.91236366, 302.91236366,
            302.91236366, 302.91236366, 302.91236366, 302.91236366);

    private List<Double> Texarray = Arrays.asList(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3);

    private List<Double> intenarray = Arrays.asList(0.482, 0.487, 0.485, 0.489, 0.485, 0.491, 0.485, 0.489, 0.489, 0.49, 0.488, 0.491, 0.487,
            0.488, 0.492, 0.487, 0.484, 0.48, 0.48, 0.478, 0.476, 0.467, 0.449, 0.42, 0.368, 0.284,
            0.187, 0.182, 0.289, 0.374, 0.411, 0.435, 0.448, 0.451, 0.458, 0.459, 0.458, 0.463, 0.456,
            0.452, 0.45, 0.44, 0.427, 0.39, 0.341, 0.236, 0.061, 0.006, 0.05, 0.246, 0.367, 0.416, 0.44,
            0.461, 0.46, 0.47, 0.475, 0.485, 0.481, 0.48, 0.483, 0.488, 0.489, 0.487, 0.492, 0.49, 0.479,
            0.476, 0.492, 0.49, 0.497, 0.5, 0.485, 0.498, 0.493, 0.492, 0.489, 0.488, 0.483, 0.477, 0.473,
            0.466, 0.449, 0.424, 0.379, 0.289, 0.173, 0.099, 0.131, 0.252, 0.329, 0.371, 0.395, 0.412, 0.42,
            0.412, 0.396, 0.394, 0.365, 0.312, 0.223, 0.086, 0.008, 0.008, 0.099, 0.255, 0.356, 0.402, 0.43,
            0.447, 0.453, 0.465, 0.477, 0.469, 0.482, 0.471, 0.481, 0.482, 0.455, 0.467, 0.466, 0.464, 0.463,
            0.461, 0.458, 0.453, 0.452, 0.438, 0.426, 0.404, 0.368, 0.305, 0.213, 0.108, 0.055, 0.06, 0.126,
            0.204, 0.247, 0.269, 0.276, 0.259, 0.23, 0.171, 0.089, 0.015, -0.008, -0.008, 0.039, 0.157, 0.257,
            0.335, 0.366, 0.399, 0.416, 0.434, 0.445, 0.449, 0.453);

    private List<Double> errarray = Arrays.asList(0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006,
            0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
            0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
            0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
            0.005, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.009, 0.007, 0.007,
            0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
            0.007, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
            0.007, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.007, 0.007, 0.007, 0.007, 0.007,
            0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
            0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
            0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
            0.003, 0.003, 0.003, 0.003, 0.003);

    private void testPars(double[] fitpars, double[] validpars, double[] fiterrs) {
        double scale = 5.0;
        for (int i = 0; i < fitpars.length; i++) {
            boolean inRange = Math.abs(fitpars[i] - validpars[i]) < fiterrs[i] * scale;
            if (!inRange) {
                System.out.println("par diff = " + i + " " + Math.abs(fitpars[i] - validpars[i]));
                System.out.println("err * " + scale + " = " + fiterrs[i] * scale);
            }
            Assert.assertTrue(inRange);
        }
    }

    private void dumpPars(String name, double[] fitpars, double[] validpars, double[] fiterrs, double[] validerrs) {
        System.out.println(name);
        System.out.print("Valid Pars:");
        for (int i = 0; i < validpars.length; i++) {
            System.out.printf("%10.3f ", validpars[i]);
        }
        System.out.print("\nFit Pars:  ");
        for (int i = 0; i < fitpars.length; i++) {
            System.out.printf("%10.3f ", fitpars[i]);
        }
        System.out.print("\nValid Err: ");
        for (int i = 0; i < validerrs.length; i++) {
            System.out.printf("%10.3f ", validerrs[i]);
        }
        System.out.print("\nFit   Err: ");
        for (int i = 0; i < fiterrs.length; i++) {
            System.out.printf("%10.3f ", fiterrs[i]);
        }
        System.out.print("\n");
    }

    private void dumpRMS(double fitRMS, double validRMS) {
        System.out.printf("Fit RMS %10.3f Valid RMS %10.3f\n", fitRMS, validRMS);
    }

    private List<Double> fixOffset(List<Double> values, double field) {
        List<Double> newValues = new ArrayList<>();
        for (double value : values) {
            newValues.add(value / (2.0 * Math.PI) / field);
        }
        return newValues;
    }

    @Test
    public void testTrottPalmer() {
        CoMDOptions options = new CoMDOptions(false);
        CESTFitter fitting = new CESTFitter(options);
        Nuclei nuc = Nuclei.findNuclei("C13");
        double field = CoMDPreferences.getRefField() * nuc.getFreqRatio();
        List<Double> fieldArray = Collections.nCopies(errarray.size(), field);
        List<Double>[] xValues = new List[4];
        xValues[0] = fixOffset(offset, field);
        xValues[1] = fixOffset(omega, 1.0);
        xValues[2] = Texarray;
        xValues[3] = fieldArray;
        fitting.setData(xValues, intenarray, errarray);

        double[] guesses = {180.795, 0.103, 2.784, -1.193, 2.383, 11.569, 88.052};

        var fitResultOpt = fitting.doFit("trott_palmer", guesses, options);
        if (fitResultOpt.isPresent()) {
            FitResult fit = fitResultOpt.get();

            double[] fitpars = fit.getPars(0);
            double validDeltaA0 = fixOffset(Arrays.asList(2705.7656176771607), field).get(0);
            double validDeltaB0 = fixOffset(Arrays.asList(-1244.26540485914), field).get(0);
            double[] validpars = new double[]{159.40549413953838, 0.07848868547481112, validDeltaA0, validDeltaB0, 2.3834525779746167,
                    2.3834525779746167, 16.071477740770437, 134.1181936680427}; //{165.02477973700172, 0.07911350221009517, 2708.726939741288, -1244.5123148064008, 2.382606761775722, 15.637464131056309, 127.51875931837206};

            double[] fiterrs = fit.getErrs(0);
            double validDeltaA0Err = fixOffset(Arrays.asList(2.1870880660250585), field).get(0);
            double validDeltaB0Err = fixOffset(Arrays.asList(4.25048585761765), field).get(0);
            double[] validerrs = new double[]{4.070771048976483, 7.910314021442713E-4, validDeltaA0Err, validDeltaB0Err, 0.004839626206279982,
                    0.004839626206279982, 0.33246477045330786, 4.989869976578456}; //{4.379924820081111, 0.0007986456868965598, 2.6001243193320365, 4.046178640951488, 0.0048026137212628375, 0.33970335899651316, 5.304491575138278};

            dumpPars("trott_palmer", fitpars, validpars, fiterrs, validerrs);

            Assert.assertArrayEquals(fiterrs, validerrs, 6.0e1);

            for (int i = 0; i < fitpars.length; i++) {
                Assert.assertTrue(Math.abs(fitpars[i] - validpars[i]) < fiterrs[i] * 4);
            }

            double fitrms = fit.getRms();
            double validrms = 0.012800369270230845;

            dumpRMS(fitrms, validrms);

            Assert.assertEquals(fitrms, validrms, 6.0e-3);
        }
    }

    @Test
    public void testSD() {
        CoMDOptions options = new CoMDOptions(false);
        CESTFitter fitting = new CESTFitter(options);
        Nuclei nuc = Nuclei.findNuclei("C13");
        double field = CoMDPreferences.getRefField() * nuc.getFreqRatio();
        List<Double> fieldArray = Collections.nCopies(errarray.size(), field);
        List<Double>[] xValues = new List[4];
        xValues[0] = fixOffset(offset, field);
        xValues[1] = fixOffset(omega, 1.0);
        xValues[2] = Texarray;
        xValues[3] = fieldArray;
        fitting.setData(xValues, intenarray, errarray);

        double[] guesses = {180.795, 0.103, 2.784, -1.193, 2.383, 11.569, 88.052};

        var fitResultOpt = fitting.doFit("sd", guesses, options);
        if (fitResultOpt.isPresent()) {
            FitResult fit = fitResultOpt.get();

            double[] fitpars = fit.getPars(0);
            double validDeltaA0 = fixOffset(Arrays.asList(2710.5848956396685), field).get(0);
            double validDeltaB0 = fixOffset(Arrays.asList(-1247.6238696374273), field).get(0);
            double[] validpars = new double[]{164.31219540872297, 0.08000559530543723, validDeltaA0, validDeltaB0, 2.3747806394904574,
                    2.3747806394904574, 16.15024270349437, 125.73252396217907};

            double[] fiterrs = fit.getErrs(0);
            double validDeltaA0Err = fixOffset(Arrays.asList(2.6506158711548755), field).get(0);
            double validDeltaB0Err = fixOffset(Arrays.asList(4.4500688239216775), field).get(0);
            double[] validerrs = new double[]{4.694876487305635, 0.0008407892782788573, validDeltaA0Err, validDeltaB0Err, 0.00510500924091006,
                    0.00510500924091006, 0.38366068154645305, 5.631002751559051};
            dumpPars("SD", fitpars, validpars, fiterrs, validerrs);

            Assert.assertArrayEquals(fiterrs, validerrs, 6.0e1);

            testPars(fitpars, validpars, fiterrs);

            double fitrms = fit.getRms();
            double validrms = 0.012583146844688359;

            dumpRMS(fitrms, validrms);

            Assert.assertEquals(fitrms, validrms, 6.0e-3);
        }
    }

    @Test
    public void testBaldwinKay() {
        CoMDOptions options = new CoMDOptions(false);
        CESTFitter fitting = new CESTFitter(options);
        Nuclei nuc = Nuclei.findNuclei("C13");
        double field = CoMDPreferences.getRefField() * nuc.getFreqRatio();
        List<Double> fieldArray = Collections.nCopies(errarray.size(), field);
        List<Double>[] xValues = new List[4];
        xValues[0] = fixOffset(offset, field);
        xValues[1] = fixOffset(omega, 1.0);
        xValues[2] = Texarray;
        xValues[3] = fieldArray;
        fitting.setData(xValues, intenarray, errarray);

        double[] guesses = {180.795, 0.103, 2.784, -1.193, 2.383, 11.569, 88.052};

        var fitResultOpt = fitting.doFit("baldwinkay", guesses, options);
        if (fitResultOpt.isPresent()) {
            FitResult fit = fitResultOpt.get();


            double[] fitpars = fit.getPars(0);
            double validDeltaA0 = fixOffset(Arrays.asList(2705.4283684664874), field).get(0);
            double validDeltaB0 = fixOffset(Arrays.asList(-1244.1504931304519), field).get(0);
            double[] validpars = new double[]{168.45940742687503, 0.07848006386286228, validDeltaA0, validDeltaB0, 2.3839104552778836,
                    2.3839104552778836, 17.777799052418256, 132.49073925170381};

            double[] fiterrs = fit.getErrs(0);
            double validDeltaA0Err = fixOffset(Arrays.asList(2.3553917781329563), field).get(0);
            double validDeltaB0Err = fixOffset(Arrays.asList(3.3808605211180436), field).get(0);
            double[] validerrs = new double[]{4.527537134628853, 0.0007171967856360292, validDeltaA0Err, validDeltaB0Err, 0.004685962985573297,
                    0.004685962985573297, 0.3932615215677069, 5.7049278404942605};

            dumpPars("BaldwinKay", fitpars, validpars, fiterrs, validerrs);

            Assert.assertArrayEquals(fiterrs, validerrs, 6.0e1);

            testPars(fitpars, validpars, fiterrs);

            double fitrms = fit.getRms();
            double validrms = 0.012936289591623533;

            dumpRMS(fitrms, validrms);

            Assert.assertEquals(fitrms, validrms, 6.0e-3);
        }
    }

    @Test
    public void testLaguerre() {
        CoMDOptions options = new CoMDOptions(false);
        CESTFitter fitting = new CESTFitter(options);
        Nuclei nuc = Nuclei.findNuclei("C13");
        double field = CoMDPreferences.getRefField() * nuc.getFreqRatio();
        List<Double> fieldArray = Collections.nCopies(errarray.size(), field);
        List<Double>[] xValues = new List[4];
        xValues[0] = fixOffset(offset, field);
        xValues[1] = fixOffset(omega, 1.0);
        xValues[2] = Texarray;
        xValues[3] = fieldArray;
        fitting.setData(xValues, intenarray, errarray);

        double[] guesses = {180.795, 0.103, 2.784, -1.193, 2.383, 88.052};

        var fitResultOpt = fitting.doFit("laguerre", guesses, options);
        if (fitResultOpt.isPresent()) {
            FitResult fit = fitResultOpt.get();

            double[] fitpars = fit.getPars(0);
            double validDeltaA0 = fixOffset(Arrays.asList(2771.0529601396183), field).get(0);
            double validDeltaB0 = fixOffset(Arrays.asList(-1285.1744405197783), field).get(0);
            double[] validpars = new double[]{288.937925578066, 0.07087265844081425, validDeltaA0, validDeltaB0, 2.4246256495623104,
                    2.4246256495623104, 8.121844525582754, 8.121844525582754};

            double[] fiterrs = fit.getErrs(0);
            double validDeltaA0Err = fixOffset(Arrays.asList(2.5674908203535574), field).get(0);
            double validDeltaB0Err = fixOffset(Arrays.asList(3.0651902549305405), field).get(0);
            double[] validerrs = new double[]{3.58822964273236, 0.0006570223547152888, validDeltaA0Err, validDeltaB0Err, 0.004205721940200411,
                    0.004205721940200411, 0.2937911953297107, 0.2937911953297107};

            dumpPars("Laguerre", fitpars, validpars, fiterrs, validerrs);

            Assert.assertArrayEquals(fiterrs, validerrs, 6.0e1);

            testPars(fitpars, validpars, fiterrs);

            double fitrms = fit.getRms();
            double validrms = 0.017326921720901785;

            dumpRMS(fitrms, validrms);

            Assert.assertEquals(fitrms, validrms, 6.0e-3);
        }
    }

//    @Test
//    public void testcestR1rhoExact1() {
//        CESTFit fitting = new CESTFit();
//        fitting.setData(offset, omega, intenarray, errarray);
//        
//        CPMGFitResult fit = fitting.doFit("cestr1rhoexact1", false, false);
//
//        double[] fitpars = fit.getPars(0);
//        double[] validpars = new double[]{169.62939271098404, 0.08068515155284373, 2712.7844655406498, -1243.9801627110282, 
//                                            2.3895199356329515, 16.787630038840575, 117.83139681295168}; 
//
//        System.out.print("cestR1rhoExact1 valid parameters: ");
//        for(int i=0; i<validpars.length; i++){
//            System.out.print(validpars[i]+" ");
//        }
//        System.out.print("\ncestR1rhoExact1 fit parameters: ");
//        for(int i=0; i<fitpars.length; i++){
//            System.out.print(fitpars[i]+" ");
//        }
//        System.out.print("\n");
//        
//        Assert.assertArrayEquals(fitpars, validpars, 6.0e-13);
//        
//        double[] fiterrs = fit.getErrs(0);
//        double[] validerrs = new double[]{4.302966390993574, 0.0008299421141931896, 2.515481529143345, 3.4487438852966332, 
//                                            0.004677381349308522, 0.35494990123302583, 4.8494380436447075}; 
//
//        System.out.print("cestR1rhoExact1 valid errors: ");
//        for(int i=0; i<validerrs.length; i++){
//            System.out.print(validerrs[i]+" ");
//        }
//        System.out.print("\ncestR1rhoExact1 fit errors: ");
//        for(int i=0; i<fiterrs.length; i++){
//            System.out.print(fiterrs[i]+" ");
//        }
//        System.out.print("\n");
//        
//        Assert.assertArrayEquals(fiterrs, validerrs, 6.0e1);
//    }
}
