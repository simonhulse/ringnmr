package org.comdnmr.modelfree;

import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.DynamicsSource;
import org.comdnmr.modelfree.models.MFModelIso;
import org.junit.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;

public class FitR1R2NOEModelTest {

    RelaxFit makeRelaxFit() {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setLambdaS(0.001);
        relaxFit.setLambdaTau(0.01);
        relaxFit.setUseLambda(false);
        relaxFit.setLogJMode(false);
        boolean fitJ = false;
        relaxFit.setFitJ(fitJ);
        return relaxFit;
    }

    MFModelIso makeModel(String modelName, double tau) {
        boolean localFitTau = false;
        double localTauFraction = 1.0;
        boolean fitExchange = false;

        MFModelIso model = MFModelIso.buildModel(modelName,
                localFitTau, tau, localTauFraction, fitExchange);
        model.setTauFraction(localTauFraction);
        return model;
    }

    R1R2NOEDataValue makeDataValue(MFModelIso model, RelaxEquations rlxEq, double[] pars, MolDataValues resData) {
        double[] valJ = model.calc(rlxEq.wValues, pars);
        double r1 = rlxEq.R1(valJ);
        double r1Error = r1 * 0.03;
        double r2 = rlxEq.R2(valJ, 0.0);
        double r2Error = r2 * 0.03;
        double noe = rlxEq.NOE(valJ);
        double noeError = 0.05;
        return new R1R2NOEDataValue(resData, r1, r1Error, r2, r2Error, noe, noeError, rlxEq);
    }

    public Map<String, MolDataValues> getMolDataValues() {
        DynamicsSource dynamicsSourceFactory = new DynamicsSource(true, true, true, true);
        double[] v = {0.0, 1.0, 2.0};
        Map<String, MolDataValues> molDataRes = new HashMap<>();
        MolDataValues resData = new MolDataValues("3.CB", v, dynamicsSourceFactory);
        molDataRes.put("a", resData);
        return molDataRes;
    }

    public PointValuePair fitModel(RelaxFit relaxFit, String modelName, double tau, double[] pars, double[] sfs) {
        MFModelIso model = makeModel(modelName, tau);

        Map<String, MolDataValues> molDataRes = getMolDataValues();
        relaxFit.setRelaxData(molDataRes);
        MolDataValues resData = molDataRes.get("a");

        for (double sf : sfs) {
            RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(sf, "H", "N");
            R1R2NOEDataValue dValue = makeDataValue(model, rlxEq, pars, resData);
            resData.addData(dValue);

        }
        double[] start = model.getStart();
        double[] lower = model.getLower();
        double[] upper = model.getUpper();
        resData.setTestModel(model);
        Optional<PointValuePair> fitResultOpt = relaxFit.fitResidueToModel(start, lower, upper);
        if (fitResultOpt.isPresent()) {
            PointValuePair fitResult = fitResultOpt.get();
            for (int i = 0; i < start.length; i++) {
                System.out.printf("%d %4.2f %4.2f %4.2f %4.2f %4.2f\n", i, start[i], lower[i], upper[i], fitResult.getPoint()[i], pars[i]);
            }
            return fitResult;
        } else {
            return null;
        }
    }

    double[] clearTau(double[] pars, double[] result) {
        if (pars[0] > 0.9999) {
            result[1] = 0.0;
        }
        if (pars[2] > 0.9999) {
            result[3] = 0.0;
        }
        return result;
    }
    @Test
    public void testModel1() {
        String modelName = "2sf";
        double tau = 10.0;
        double[] pars = {0.3, 0.0, 0.3, 0.2};
        double[] sfs = {600.0e6, 800.0e6};
        RelaxFit relaxFit = makeRelaxFit();

        PointValuePair fitResult = fitModel(relaxFit, modelName, tau, pars, sfs);
        Assert.assertArrayEquals(pars, fitResult.getPoint(), 0.01);
    }

    @Test
    public void testModel2() {
        String modelName = "2sf";
        double tau = 100.0;
        double[] pars = {0.3, 0.0, 0.3, 0.2};
        double[] sfs = {600.0e6, 800.0e6};

        RelaxFit relaxFit = makeRelaxFit();
        PointValuePair fitResult = fitModel(relaxFit, modelName, tau, pars, sfs);
        Assert.assertArrayEquals(pars, fitResult.getPoint(), 0.01);
    }

    @Test
    public void testModel3() {
        String modelName = "2sf";
        double tau = 10.0;
        double[] pars = {0.9, 0.0, 0.3, 0.2};
        double[] sfs = {600.0e6, 800.0e6};

        RelaxFit relaxFit = makeRelaxFit();
        PointValuePair fitResult = fitModel(relaxFit, modelName, tau, pars, sfs);
        Assert.assertArrayEquals(pars, fitResult.getPoint(), 0.01);
    }
    @Test
    public void testModel4() {
        String modelName = "2sf";
        double tau = 10.0;
        double[] pars = {0.9, 0.1, 1.0, 0.0};
        double[] sfs = {600.0e6, 800.0e6};

        RelaxFit relaxFit = makeRelaxFit();
        PointValuePair fitResult = fitModel(relaxFit, modelName, tau, pars, sfs);
        double[] resultJ = clearTau(pars, fitResult.getPoint());
        Assert.assertArrayEquals(pars,resultJ, 0.01);
    }

    public void testJCalc() {
        String modelName = "2sf";
        double tau = 100.0;
        double[] pars = {0.3, 0.0, 0.3, 0.2};
        double sf = 500.0e6;
        Map<String, MolDataValues> molDataRes = getMolDataValues();
        MolDataValues resData = molDataRes.get("a");

        MFModelIso model = makeModel(modelName, tau);
        RelaxEquations rlxEq = RelaxEquations.getRelaxEquations(sf, "H", "N");
        R1R2NOEDataValue dValue = makeDataValue(model, rlxEq, pars, resData);
        double[] valJ = model.calc(rlxEq.wValues, pars);

        SpectralDensityCalculator spectralDensityCalculator = new SpectralDensityCalculator();

        PointValuePair result = spectralDensityCalculator.fit(dValue);
        System.out.println("sd fit " + result);
        double[] jResult = result.getPoint();
        for (int k = 0; k < jResult.length; k++) {
            System.out.printf("%.3g %.3g\n", valJ[k], jResult[k]);
        }
        Assert.assertArrayEquals(valJ, result.getPoint(), 1.0e-12);
    }
}