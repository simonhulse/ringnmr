package org.comdnmr.modelfree.models;

import org.comdnmr.modelfree.RelaxEquations;
import org.junit.Assert;
import org.junit.Test;

public class MFModelTest {
    public double[][] getJForModels(String modelName, double[] pars2sf, double[] parsTest) {
        double sf = 500.0e8;
        boolean localFitTau =false;
        double tau = 100.0;
        double localTauFraction = 1.0;
        boolean fitExchange = false;


        MFModelIso model2sf = MFModelIso.buildModel("2sf",
                localFitTau, tau, localTauFraction, fitExchange);

        RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(sf, "H", "N");
        double[] valJ2sf = model2sf.calc(rlxEq.getW(), pars2sf);

        MFModelIso model1f = MFModelIso.buildModel(modelName,
                localFitTau, tau, localTauFraction, fitExchange);
        double[] valJ1f = model1f.calc(rlxEq.getW(), parsTest);
        return new double[][]{valJ2sf, valJ1f};

    }
    @Test
    public void test1f() {
        double[] parsTest = {0.8, 0.1};
        double[] pars2sf = {0.8, 0.1, 1.0, 0.0};
        double[][] jValues = getJForModels("1f", pars2sf, parsTest);
        for (int i = 0;i< jValues[0].length;i++) {
            double delta = Math.abs(jValues[0][i] * 1.0e-3);
            System.out.println(jValues[0][i] + " " + delta);
            Assert.assertEquals(jValues[0][i], jValues[1][i], delta);
        }
    }
    @Test
    public void test1s() {
        double[] parsTest = {0.8, 1.0};
        double[] pars2sf = {0.8, 0.1, 1.0, 0.0};
        double[][] jValues = getJForModels("1s", pars2sf, parsTest);
        Assert.assertArrayEquals(jValues[0], jValues[1], 1.0e-9);
    }
    @Test
    public void test2s() {
        double[] parsTest = {0.7, 0.9, 1.0};
        double[] pars2sf = {0.7, 0.0, 1.0, 0.9};
        double[][] jValues = getJForModels("2s", pars2sf, parsTest);
        Assert.assertArrayEquals(jValues[0], jValues[1], 1.0e-9);
    }

}