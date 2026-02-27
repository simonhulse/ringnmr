/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.DataIO;
import org.comdnmr.data.DynamicsSource;
import org.comdnmr.modelfree.*;
import org.comdnmr.modelfree.RelaxFit.DiffusionType;
import org.comdnmr.modelfree.models.*;
import org.comdnmr.util.CoMDPreferences;
import org.junit.Test;
import org.junit.Assert;
import org.nmrfx.chemistry.MoleculeFactory;

/**
 *
 * @author Martha
 */
public class DRefineTest {

//    private final double[] guesses = {0.75 * Di, Di, 1.25 * Di, Math.PI / 2, Math.PI / 2, Math.PI / 2};
//    private final double[] guesses = {4.4170 * 1e7, 4.5832 * 1e7, 6.0129 * 1e7, Math.toRadians(98.06),
//            Math.toRadians(68.64), Math.toRadians(77.42)};
    private final double[] fields = {400.0e6, 500e6, 600e6, 700e6, 800.0e6};

    private Map<String, MolDataValues> loadVectors(String vectorFileName) throws IOException {
        File file = new File(vectorFileName);
        Path path = file.toPath();
        Map<String, MolDataValues> molValues = new TreeMap<>();
        DynamicsSource dynamicsSourceFactory = new DynamicsSource(true, true, true, true);
        CoMDPreferences.setOptimizer("CMAES");
        Files.lines(path).forEach(line -> {
            String[] fields = line.split("\t");
            if (fields.length == 5) {
                int iRes = Integer.valueOf(fields[0]);
                String chain = "A";
                String atomSpecifier = chain + ":" + iRes + ".H";
                double[] v = new double[3];
                v[0] = Double.valueOf(fields[2]);
                v[1] = Double.valueOf(fields[3]);
                v[2] = Double.valueOf(fields[4]);
                MolDataValues molDataValue = new MolDataValues(atomSpecifier, v,
                        dynamicsSourceFactory);
                molValues.put(atomSpecifier, molDataValue);
            }
        });

        return molValues;
    }

    private void loadData(String relaxDataFileName, Map<String, MolDataValues> molDataValues) throws IOException {

        File relaxFile = new File(relaxDataFileName);
        Path relaxPath = relaxFile.toPath();
        Map<Integer, Integer> fieldMap = new HashMap<>();
        Files.lines(relaxPath).forEach(line -> {
            String[] fields = line.split("\\s+");
            if (fields.length == 11) {
                int iRes = Integer.valueOf(fields[0]);
                String chain = fields[1];
                String resSpecifier = chain + ":" + iRes + ".H";
                MolDataValues molData = molDataValues.get(resSpecifier);
                if (molData == null) {
                    System.out.println("no res " + resSpecifier);
                } else {
                    double field = Double.valueOf(fields[4]);
                    String e1 = fields[3];
                    String e2 = fields[2];
                    RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, e1, e2);

                    double r1 = Double.valueOf(fields[5]);
                    double r1Error = Double.valueOf(fields[6]);
                    double r2 = Double.valueOf(fields[7]);
                    double r2Error = Double.valueOf(fields[8]);
                    double noe = Double.valueOf(fields[9]);
                    double noeError = Double.valueOf(fields[10]);
                    R1R2NOEDataValue dValue = new R1R2NOEDataValue(molData, r1, r1Error, r2, r2Error, noe, noeError, relaxObj);
                    molData.addData(dValue);
                }
            }
        });
        Iterator<Map.Entry<String, MolDataValues>> iter = molDataValues.entrySet().iterator();
        while (iter.hasNext()) {
            Map.Entry<String, MolDataValues> entry = iter.next();
            if (entry.getValue().getData().isEmpty()) {
                iter.remove();
            }

        }
    }

    public double[] scalePars(double[] pars, int nDiffPars) {
        double[] scaledPars = new double[pars.length];
        for (int p = 0; p < pars.length; p++) {
            double newPar = pars[p];
            if (p < nDiffPars) {
                newPar /= 1.0e7;
            } else if (p >= nDiffPars) {
                newPar *= 180.0 / Math.PI;
            }
            scaledPars[p] = newPar;
        }
        return scaledPars;
    }

    private Map<String, MolDataValues> loadTestData() {
        try {
            Map<String, MolDataValues> molDataValues = loadVectors("src/test/data/1P7F_A_vectors.txt");
            loadData("src/test/data/1P7F_data.txt", molDataValues);
            return molDataValues;
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
            return null;
        }

    }

    private Map<String, MolDataValues> loadDeuteriumTestData() {
        File file = new File("src/test/data/val121gamma2no700.csv");
        try {
            DataIO.loadRelaxationTextFile(file);
        } catch (IOException e) {
            return null;
        }
        return FitDeuteriumModel.getData(false);
    }

    @Test
    public void testDeuteriumD1() {
        MoleculeFactory.clearAllMolecules();
        MoleculeFactory.setActive(null);
        var data = loadDeuteriumTestData();
        FitDeuteriumModel fitDeuteriumModel = new FitDeuteriumModel();
        fitDeuteriumModel.setTau(10.0);
        var modelNames = List.of("1f");

        fitDeuteriumModel.testModels(data,modelNames);
    }

    @Test
    public void testDeuteriumD1J() {
        MFModelIso1 model = new MFModelIso1();
        model.setSScale(9.0);

        double[] hFreqs = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0};
        double[] omegas = new double[hFreqs.length * 2 + 1];

        for (int i = 0; i < hFreqs.length; i++) {
            omegas[i * 2 + 1] = 1.0e6* hFreqs[i] * RelaxEquations.GAMMA_D / RelaxEquations.GAMMA_H * 2.0 * Math.PI;
            omegas[i * 2 + 2] = omegas[i * 2 + 1] * 2;
        }
        Arrays.sort(omegas);
        double[] pars = {7.0, 1.0};
        var jValues = model.calc(omegas, pars);
        for (int i = 0; i < omegas.length; i++) {
            System.out.printf("%9.4f %9.4f\n", omegas[i]*1.0e-9, Math.log10(jValues[i]*1.0e9));
        }
    }

    @Test
    public void testDeuteriumD1fJ() {
        MFModelIso1f model = new MFModelIso1f();
        model.setSScale(9.0);
        double[] hFreqs = {100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0};
        double[] omegas = new double[hFreqs.length * 2 + 1];

        for (int i = 0; i < hFreqs.length; i++) {
            omegas[i * 2 + 1] = 1.0e6* hFreqs[i] * RelaxEquations.GAMMA_D / RelaxEquations.GAMMA_H * 2.0 * Math.PI;
            omegas[i * 2 + 2] = omegas[i * 2 + 1] * 2;
        }
        Arrays.sort(omegas);
        double[] pars = {12.411, 0.471, 44.9e-3};
        var jValues = model.calc(omegas, pars);
        for (int i = 0; i < omegas.length; i++) {
            System.out.printf("%9.4f %9.4f\n", omegas[i]*1.0e-9, Math.log10(jValues[i]*1.0e9));
        }
    }

    public void testIsoModel(RelaxFit relaxFit, int modelNum) {
        DRefineTest.this.testIsoModel(relaxFit, modelNum, null);
    }

    public void testIsoModel(RelaxFit relaxFit, int modelNum, String matchSpec) {
        double tau = 5.0e-9;
        MFModelIso model;
//            public MFModelIso1(boolean fitTau, double targetTau, double tauFraction,
//            boolean includeEx) {

        switch (modelNum) {
            case 1:
                model = new MFModelIso1(!relaxFit.isUseGlobalTau(), tau, 0.1, false);
                break;
            case 2:
                model = new MFModelIso1f(!relaxFit.isUseGlobalTau(), tau, 0.1, false);
                break;
            case 5:
                model = new MFModelIso2s(!relaxFit.isUseGlobalTau(), tau, 0.1, false);
                break;
            case 6:
                model = new MFModelIso2sf(!relaxFit.isUseGlobalTau(), tau, 0.1, false);
                break;
            default:
                model = new MFModelIso1(!relaxFit.isUseGlobalTau(), tau, 0.1, false);
        }
        DRefineTest.this.testIsoModel(relaxFit, model, matchSpec);
    }

    public void testIsoModel(RelaxFit relaxFit, MFModelIso model, String matchSpec) {
        Map<String, MolDataValues> molData = loadTestData();
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        for (String key : molData.keySet()) {
            if ((matchSpec != null) && !key.equals(matchSpec)) {
                continue;
            }
            molDataRes.clear();
            MolDataValues resData = molData.get(key);
            if (!resData.getData().isEmpty()) {
                resData.setTestModel(model);
                molDataRes.put(key, molData.get(key));
                relaxFit.setRelaxData(molDataRes);
                double[] start = model.getStart();
                double[] lower = model.getLower();
                double[] upper = model.getUpper();
                Optional<PointValuePair> fitOpt = relaxFit.fitResidueToModel(start, lower, upper);
                fitOpt.ifPresent(fitResult -> {
                    double[] values = fitResult.getPoint();
                    double score = fitResult.getValue();
                    for (double val : values) {
                        System.out.print(val + " ");
                    }
                    System.out.println(score + " " + key);
                });
            }
        }
    }

    @Test
    public void testModel6OneRes() {
        System.out.println("test constrain 5 res");

        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setUseGlobalTau(true);
        relaxFit.setGlobalTau(3.27e-9);
        DRefineTest.this.testIsoModel(relaxFit, 2, "A44");
    }

    @Test
    public void testModel1() {
        System.out.println("test 1");
        RelaxFit relaxFit = new RelaxFit();
        testIsoModel(relaxFit, 1);
    }

    @Test
    public void testModel6Constrain() {
        System.out.println("test constrain 6");

        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setUseGlobalTau(true);
        relaxFit.setGlobalTau(3.27e-9);
        testIsoModel(relaxFit, 6);
    }

    @Test
    public void testModel2() {
        System.out.println("test 2");
        RelaxFit relaxFit = new RelaxFit();
        testIsoModel(relaxFit, 2);
    }

    @Test
    public void testModel5() {
        System.out.println("test 5");
        RelaxFit relaxFit = new RelaxFit();
        testIsoModel(relaxFit, 5);
    }

    @Test
    public void testModel6() {
        System.out.println("test 6");
        RelaxFit relaxFit = new RelaxFit();
        testIsoModel(relaxFit, 6);
    }

    @Test
    public void testMultiResidueModel() {
        System.out.println("test multi 1");

        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setUseGlobalTau(true);
        Map<String, MolDataValues> molData = loadTestData();
        double tau = 5.0e-9;
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        int nModelPars = 0;
        for (String key : molData.keySet()) {
            MolDataValues resData = molData.get(key);
            MFModelIso model = new MFModelIso1(false, 0.1, 0.0, false);
            if (!resData.getData().isEmpty()) {
                nModelPars += model.getNPars();
                resData.setTestModel(model);
                molDataRes.put(key, molData.get(key));
            }
        }
        int nGuesses = 1 + nModelPars;

        double[] guesses = new double[nGuesses];
        double[] lower = new double[nGuesses];
        double[] upper = new double[nGuesses];

        int start = 1;
        for (String key : molDataRes.keySet()) {
            MolDataValues resData = molData.get(key);
            MFModelIso model = (MFModelIso) resData.getTestModel();
            double[] resStart = model.getStart();
            double[] resLower = model.getLower();
            double[] resUpper = model.getUpper();

            System.arraycopy(resStart, 0, guesses, start, resStart.length);
            System.arraycopy(resLower, 0, lower, start, resStart.length);
            System.arraycopy(resUpper, 0, upper, start, resStart.length);
            start += resStart.length;
        }
        guesses[0] = tau;
        lower[0] = tau / 4.0;
        upper[0] = tau * 4.0;
        relaxFit.setRelaxData(molDataRes);
        PointValuePair fitResult = relaxFit.fitMultiResidueToModel(guesses, lower, upper);
        double[] values = fitResult.getPoint();
        double score = fitResult.getValue();
        for (double val : values) {
            System.out.println(val + " ");
        }
        System.out.println(score);

    }

    @Test
    public void testValueDMat1() {
        RelaxFit relaxFit = new RelaxFit();
        Map<String, MolDataValues> molData = loadTestData();
        for (MolDataValues molDataValue : molData.values()) {
            MFModelAniso1 model = new MFModelAniso1(DiffusionType.ANISOTROPIC, molDataValue.getVector());
            molDataValue.setTestModel(model);
        }
        relaxFit.setRelaxData(molData);
        relaxFit.setDiffusionType(DiffusionType.ANISOTROPIC);
        double[] pars = {4.4170 * 1e7, 4.5832 * 1e7, 6.0129 * 1e7, Math.toRadians(98.06),
            Math.toRadians(68.64), Math.toRadians(77.42)};

        double value = relaxFit.valueDMat(pars, null);
        System.out.println("value " + value);
        Assert.assertEquals(0.0, value, 0.11);
    }

    @Test
    public void testValueDMatFile() {
        RelaxFit relaxFit = new RelaxFit();
        Map<String, MolDataValues> molData = loadTestData();
        relaxFit.setRelaxData(molData);
        double tauCGuess = 3.3e-9;

        double[] rotDifPars = {4.4170, 4.5832, 6.0129, 98.06, 68.64, 77.42};
        double[] deltas = {0.125, 0.125, 0.125, 2.0, 6.0, 50.0};
        double bestFitRMS = Double.MAX_VALUE;
        double[] bestGuesses = null;
        double[] bestFitPars = null;
        DiffusionType bestType = null;
        double isoD = 1.0 / (6.0 * tauCGuess);
        for (DiffusionType diffType : DiffusionType.values()) {
            for (MolDataValues molDataValue : molData.values()) {
                MFModel model;
                if (diffType == DiffusionType.ISOTROPIC) {
                    model = new MFModelIso1();
                } else {
                    model = new MFModelAniso1(diffType, molDataValue.getVector());
                }
                molDataValue.setTestModel(model);
            }
            int nPars = diffType.getNAnglePars() + diffType.getNDiffusionPars();
            double[] guess = new double[nPars];
            relaxFit.setDiffusionType(diffType);
            System.arraycopy(diffType.getGuess(isoD), 0, guess, 0, diffType.getNDiffusionPars());
            int nGuesses = diffType.getNAngleGuesses();
            for (int iGuess = 0; iGuess < nGuesses; iGuess++) {
                double[] angleGuess = diffType.getAngles(iGuess);
                System.arraycopy(angleGuess, 0, guess, diffType.getNDiffusionPars(), diffType.getNAnglePars());
                PointValuePair fitResult = relaxFit.fitDiffusion(guess);
                double fitRMS = fitResult.getValue();
                double[] scaledGuesses = scalePars(guess, diffType.getNDiffusionPars());
                double[] fitPars = relaxFit.getPars();
                Arrays.sort(fitPars, 0, diffType.getNDiffusionPars());
                double[] scaledPars = scalePars(fitPars, diffType.getNDiffusionPars());

                System.out.println("RMS: " + fitRMS + " " + diffType + " guess " + iGuess);
                System.out.println("Scaled guesses, Fit Pars, RotDif Pars: ");
                for (int i = 0; i < scaledGuesses.length; i++) {
                    System.out.printf("guess %7.3f best Fit %7.3f RotDif %7.3f\n",
                            scaledGuesses[i], scaledPars[i], rotDifPars[i]);
                }
                System.out.println();
                if (fitRMS < bestFitRMS) {
                    bestGuesses = guess.clone();
                    bestFitRMS = fitRMS;
                    bestFitPars = fitPars.clone();
                    bestType = diffType;
                }
            }
        }
        relaxFit.setDiffusionType(bestType);
        double[][] D = relaxFit.parsToD(bestFitPars, bestType);
        double[][] VT = relaxFit.parsToVT(bestFitPars, bestType);

        double[] scaledBestGuesses = scalePars(bestGuesses, bestType.getNDiffusionPars());
        double[] scaledBestFitPars = scalePars(bestFitPars, bestType.getNDiffusionPars());

        System.out.println("\nbest Fit RMS: " + bestFitRMS + " " + bestType.toString());
        System.out.println("Scaled best guesses, best Fit Pars, RotDif Pars: ");
        for (int i = 0; i < scaledBestGuesses.length; i++) {
            System.out.printf("guess %7.3f best Fit %7.3f RotDif %7.3f\n",
                    scaledBestGuesses[i], scaledBestFitPars[i], rotDifPars[i]);
        }
        System.out.println();
        Assert.assertEquals(0.1, bestFitRMS, 0.01);
        for (int i = 0; i < scaledBestGuesses.length; i++) {
            Assert.assertEquals(rotDifPars[i], scaledBestFitPars[i], deltas[i]);
        }
        Assert.assertEquals(bestType, DiffusionType.ANISOTROPIC);
        //  relaxFit.dumpValues(bestFitPars);

    }

}
