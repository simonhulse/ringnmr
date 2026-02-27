package org.comdnmr.modelfree;

import java.util.*;

import org.apache.commons.math3.geometry.euclidean.threed.NotARotationMatrixException;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationConvention;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.Fitter;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.ANISOTROPIC;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.OBLATE;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.PROLATE;
import org.comdnmr.modelfree.models.MFModel;
import org.comdnmr.modelfree.models.MFModelAniso;
import org.comdnmr.modelfree.models.MFModelAniso1;
import org.comdnmr.modelfree.models.MFModelAniso2;
import org.comdnmr.modelfree.models.MFModelAniso5;
import org.comdnmr.modelfree.models.MFModelAniso6;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso1;
import org.comdnmr.modelfree.models.MFModelIso1f;
import org.comdnmr.modelfree.models.MFModelIso2s;
import org.comdnmr.modelfree.models.MFModelIso2sf;

/**
 *
 * @author brucejohnson
 */
public class RelaxFit {

    double lambdaS = 0.0;
    double lambdaTau = 0.0;
    boolean useLambda = false;
    boolean logJMode = false;
    boolean fitJ = false;
    Map<String, MolDataValues> molDataValues;
    double[] bestPars;
    double[] parErrs;
    double bestAIC;
    double bestChiSq;
    Fitter bestFitter;
    int[] nParsPerModel = {0, 1, 2, 0, 0, 3, 4};
    DiffusionType diffusionType;
    static final int[][] ANISO_ANGLE_STARTS = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1}
    };
    static final int[][] ANGLE_STARTS = {
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1}
    };

    double globalTau = 4.0e-9;
    boolean useGlobalTau = false;

    public double getGlobalTau() {
        return globalTau;
    }

    public void setGlobalTau(double globalTau) {
        this.globalTau = globalTau;
    }

    public boolean isUseGlobalTau() {
        return useGlobalTau;
    }

    public void setUseGlobalTau(boolean useGlobalTau) {
        this.useGlobalTau = useGlobalTau;
    }

    public boolean getFitJ() {
        return fitJ;
    }

    public void setFitJ(boolean value) {
        fitJ = value;
    }
    public void setLogJMode(boolean value) {
        logJMode = value;
    }

    public boolean getLogJMode() {
        return logJMode;
    }

    public double getLambdaS() {
        return useLambda ? lambdaS : 0.0;
    }

    public void setLambdaS(double value) {
        this.lambdaS = value;
    }

    public double getLambdaTau() {
        return useLambda ? lambdaTau : 0.0;
    }

    public void setLambdaTau(double value) {
        this.lambdaTau = value;
    }

    public void setUseLambda(boolean value) {
        this.useLambda = value;
    }

    public boolean useLambda() {
        return useLambda && (lambdaS > 1.0e-8 || lambdaTau > 1.0e-8);
    }

    public static double[] getDValues(double isoD) {
        return new double[]{0.75 * isoD, isoD, 1.25 * isoD};
    }

    public void setRelaxData(Map<String, MolDataValues> molDataValues) {
        this.molDataValues = molDataValues;
    }

    public void setDiffusionType(DiffusionType type) {
        this.diffusionType = type;
    }

    public enum DiffusionType {
        ISOTROPIC(1, 0) {
            @Override
            public double[] getGuess(double isoD) {
                return new double[]{isoD};
            }

            public double[] getAngles(int iStart) {
                return new double[0];
            }
        },
        PROLATE(2, 2) {
            @Override
            public double[] getGuess(double isoD) {
                return new double[]{0.75 * isoD, 1.25 * isoD};
            }

            public double[] getAngles(int iStart) {
                int[] jAng = ANGLE_STARTS[iStart];
                double[] angles = new double[2];
                for (int i = 0; i < 2; i++) {
                    angles[i] = (Math.PI * (jAng[i] * 2 + 1)) / 4;
                }
                return angles;
            }
        }, OBLATE(2, 2) {
            @Override
            public double[] getGuess(double isoD) {
                return new double[]{0.75 * isoD, 1.25 * isoD};
            }

            public double[] getAngles(int iStart) {
                int[] jAng = ANGLE_STARTS[iStart];
                double[] angles = new double[2];
                for (int i = 0; i < 2; i++) {
                    angles[i] = (Math.PI * (jAng[i] * 2 + 1)) / 4;
                }
                return angles;
            }
        }, ANISOTROPIC(3, 3) {
            @Override
            public double[] getGuess(double isoD) {
                return new double[]{0.75 * isoD, isoD, 1.25 * isoD};
            }

            public double[] getAngles(int iStart) {
                int[] jAng = ANISO_ANGLE_STARTS[iStart];
                double[] angles = new double[3];
                for (int i = 0; i < 3; i++) {
                    angles[i] = (Math.PI * (jAng[i] * 2 + 1)) / 4;
                }
                return angles;
            }
        };

        final int nAnglePars;
        final int nDiffPars;
        final int nAngleGuesses;

        DiffusionType(int nDiffPars, int nAnglePars) {
            this.nAnglePars = nAnglePars;
            this.nDiffPars = nDiffPars;
            this.nAngleGuesses = (int) Math.pow(2, nAnglePars);
        }

        public abstract double[] getGuess(double isoD);

        public abstract double[] getAngles(int iStart);

        public int getNAnglePars() {
            return nAnglePars;
        }

        public int getNDiffusionPars() {
            return nDiffPars;
        }

        public int getNAngleGuesses() {
            return nAngleGuesses;
        }

    }

    public double[] getJ(double[] pars, RelaxEquations relaxObj, int modelNum) {
        double tauM = pars[0];//4.5e-9;
        double s2 = pars[1];
        double[] valJ = new double[5];
        switch (modelNum) {
            case 1:
                MFModelIso1 model1 = new MFModelIso1(tauM);
                valJ = model1.calc(relaxObj.wValues, s2);
                break;
            case 2:
                double tau = pars[2];
                MFModelIso1f model2 = new MFModelIso1f(tauM);
                valJ = model2.calc(relaxObj.wValues, s2, tau);
                break;
            case 5:
                tau = pars[2];
                double sf2 = pars[3];
                MFModelIso2s model5 = new MFModelIso2s(tauM);
                valJ = model5.calc(relaxObj.wValues, s2, tau, sf2);
                break;
            case 6:
                tau = pars[2];
                sf2 = pars[3];
                double tauS = pars[4];
                MFModelIso2sf model6 = new MFModelIso2sf(tauM);
                valJ = model6.calc(relaxObj.wValues, s2, tau, sf2, tauS);
                break;
            default:
                break;
        }
        return valJ;
    }

    public double[][] parsToD(double[] pars, DiffusionType dType) {
        double Dxx = pars[0];
        double Dyy;
        double Dzz;
        switch (dType) {
            case PROLATE:
                //Dxx = Dyy
                Dyy = pars[0];
                Dzz = pars[1];
                break;
            case OBLATE:
                //Dyy = Dzz
                Dyy = pars[1];
                Dzz = pars[1];
                break;
            case ANISOTROPIC:
                Dyy = pars[1];
                Dzz = pars[2];
                break;
            default:
                Dyy = pars[0];
                Dzz = pars[0];

        }
        return new double[][]{{Dxx, 0.0, 0.0},
        {0.0, Dyy, 0.0},
        {0.0, 0.0, Dzz}};
    }

    public double[][] parsToVT(double[] pars, DiffusionType dType) {
        Rotation rot = getDRotation(pars, dType);
        return getRotationMatrix(rot);
    }

    public double[] getJDiffusion(double[] pars, RelaxEquations relaxObj, MFModelAniso model, double[] v, double[][] D, double[][] VT) {
        int extraParStart = diffusionType.getNDiffusionPars() + diffusionType.getNAnglePars();
        double[] valJ;
        double[] modelPars = new double[model.getNPars()];
        System.arraycopy(pars, extraParStart, modelPars, 0, model.getNPars());
        model.update(D, VT);
        valJ = model.calc(relaxObj.wValues, modelPars);
        return valJ;
    }

    public double[] getJDiffusion(double[] pars, RelaxEquations relaxObj, int modelNum, double[] v, double[][] valD, double[][] valVT) {
        int extraParStart = diffusionType.getNDiffusionPars() + diffusionType.getNAnglePars();
        double s2 = pars[extraParStart];
        double[] valJ = new double[5];
        switch (modelNum) {
            case 1:
                MFModelAniso1 model1 = new MFModelAniso1(diffusionType, valD, valVT, v);
                valJ = model1.calc(relaxObj.wValues, s2);
                break;
            case 2:
                double tau = pars[extraParStart + 1];
                MFModelAniso2 model2 = new MFModelAniso2(diffusionType, valD, valVT, v);
                valJ = model2.calc(relaxObj.wValues, s2, tau);
                break;
            case 5:
                tau = pars[extraParStart + 1];
                double sf2 = pars[extraParStart + 2];
                MFModelAniso5 model5 = new MFModelAniso5(diffusionType, valD, valVT, v);
                valJ = model5.calc(relaxObj.wValues, s2, tau, sf2);
                break;
            case 6:
                tau = pars[extraParStart + 1];
                sf2 = pars[extraParStart + 2];
                double tauS = pars[extraParStart + 3];
                MFModelAniso6 model6 = new MFModelAniso6(diffusionType, valD, valVT, v);
                valJ = model6.calc(relaxObj.wValues, s2, tau, sf2, tauS);
                break;
            default:
                break;
        }
        return valJ;
    }

    public Rotation getDRotation(double[] pars, DiffusionType dType) {
        int nEqlDiffPars;
        double alpha;
        double beta;
        double gamma;
        switch (dType) {
            case PROLATE, OBLATE -> {
                //Dyy = Dzz
                //Dxx = Dyy
                nEqlDiffPars = 1;
                alpha = pars[3 - nEqlDiffPars];
                beta = pars[4 - nEqlDiffPars];
                gamma = 0;
            }
            default -> {
                alpha = pars[3];
                beta = pars[4];
                gamma = pars[5];
            }
        }
//        System.out.println("alpha = " + alpha*180./Math.PI + " beta = " + beta*180./Math.PI + " gamma = " + gamma*180./Math.PI);
        Rotation rot = null;
        try {
            rot = new Rotation(RotationOrder.ZYZ, RotationConvention.VECTOR_OPERATOR, alpha, beta, gamma);
        } catch (NotARotationMatrixException nE) {
            System.out.println("Can't create rot mat:" + nE.getMessage());
            double[][] rotMatCatch = rot.getMatrix();
            for (int i = 0; i < 3; i++) {
                rotMatCatch[1][i] = -rotMatCatch[1][i];
            }
            try {
                rot = new Rotation(rotMatCatch, 1e-6);
            } catch (NotARotationMatrixException nE2) {
                System.out.println("Can't create rot mat 2nd try:" + nE.getMessage());
                rot = null;
            }
        }
        return rot;
    }

    public double[][] getRotationMatrix(Rotation rot) {
        return new Array2DRowRealMatrix(rot.getMatrix()).transpose().getData();
    }

    double[] calcDeltaSqJ(MolDataValues molData, double[] resPars, MFModel testModel, boolean report) {
        double sumSq = 0.0;
        double sumSqNW = 0.0;
        double[][] jValues = molData.getJValues();
        double[] jCalc = testModel.calc(jValues[0], resPars);
        int nValues = jCalc.length;
        double[] weights = jValues[jValues.length - 1];
        for (int i=0;i< jCalc.length;i++) {
            final double delta;
            final double delta2;
            final double jErr;
            if (logJMode) {
                delta = Math.log10(jCalc[i]) - Math.log10(jValues[1][i]);
                delta2 = delta * delta;
                double high = jValues[1][i] + jValues[2][i];
                double low = jValues[1][i] - jValues[2][i];
                jErr = Math.abs(Math.log10(high) - Math.log10(low)) / 2.0;
                double jErr2 = jErr * jErr;
                sumSq += weights[i] * delta2 / jErr2;
                sumSqNW += delta2;
            } else {
                delta = jCalc[i] * 1.0e9 - jValues[1][i] * 1.0e9;
                delta2 = delta * delta;
                jErr = jValues[2][i] * 1.0e9;
                double jErr2 = jErr * jErr;
                sumSq += weights[i] * delta2 / jErr2;
                sumSqNW += delta2;
            }

         if (report) {
             System.out.printf("%3d %11.5g %11.5g %11.5g %11.5g %11.5g %11.5g\n", i, jValues[0][i] * 1.0e-9, jCalc[i] * 1.0e9, jValues[1][i] * 1.0e9, jErr, weights[i], delta);
         }
        }
       if (report) {
           System.out.printf("%11.5g %11.5g\n",Math.sqrt(sumSq/jCalc.length), Math.sqrt(sumSqNW/jCalc.length));
       }
        double complexityS = testModel.getComplexityS();
        double complexityTau = testModel.getComplexityTau();
        return new double[]{sumSq, complexityS, complexityTau, nValues};
    }

    double[] calcDeltaSqR(MolDataValues molData, double[] resPars, MFModel testModel) {
        double sumComplexityS = 0.0;
        double sumComplexityTau = 0.0;
        double sumSq = 0.0;
        int nPar = 0;
        for (RelaxDataValue value : molData.getData()) {
            R1R2NOEDataValue dValue  = (R1R2NOEDataValue) value;
            RelaxEquations relaxObj = dValue.relaxObj;
            double[] J = testModel.calc(relaxObj.wValues, resPars);
            sumComplexityS += testModel.getComplexityS();
            sumComplexityTau += testModel.getComplexityTau();
            double r1 = relaxObj.R1(J);
            // fixme rEx should be field dependent
            double rEx = testModel.includesEx() ? resPars[resPars.length - 1] : 0.0;
            double r2 = relaxObj.R2(J, rEx);
            double noe = relaxObj.NOE(J);
            double delta2 = dValue.score2(r1, r2, noe);
            sumSq += delta2;
            nPar += 3;
        }
        return new double[]{sumSq, sumComplexityS, sumComplexityTau, nPar};
    }

    double[] calcDeltaSq(MolDataValues molData, double[] resPars, MFModel testModel, boolean report) {
        if (fitJ) {
            return calcDeltaSqJ(molData, resPars, testModel, report);
        } else {
            return calcDeltaSqR(molData, resPars, testModel);
        }
    }

    public Score score(double[] pars, boolean keepPars) {
        return score(pars, keepPars, false);
    }

    public Score score(double[] pars, boolean keepPars, boolean report) {
        double sumSq = 0.0;
        int n = 0;
        int nComplex = 0;
        boolean parsOK = true;
        double sumComplexityS = 0.0;
        double sumComplexityTau = 0.0;
        for (MolDataValues molData : molDataValues.values()) {
            MFModel testModel = molData.getTestModel();
            double[] resPars;
            if (useGlobalTau) {
                int nResPars = testModel.getNPars();
                resPars = new double[nResPars + 1];
                resPars[0] = globalTau;
                System.arraycopy(pars, 0, resPars, 1, nResPars);
            } else {
                resPars = pars;
            }
            double[] resResult = calcDeltaSq(molData, resPars, testModel, report);
            sumSq += resResult[0];
            sumComplexityS += resResult[1];
            sumComplexityTau += resResult[2];
            if (fitJ) {
                nComplex++;
            } else {
                nComplex += molData.getData().size();
            }
            n += (int) Math.round(resResult[3]);

            if (!testModel.checkParConstraints()) {
                parsOK = false;
            }
        }
        double avgComplexityS = sumComplexityS / nComplex;
        double avgComplexityTau = sumComplexityTau / nComplex;
        Score score;
        if (keepPars) {
            score = new Score(sumSq, n, pars.length, parsOK, avgComplexityS, avgComplexityTau, pars.clone());
        } else {
            score = new Score(sumSq, n, pars.length, parsOK, avgComplexityS, avgComplexityTau);
        }
        return score;
    }

    public Map<String, MolDataValues> genBootstrap(Random random, MFModel model, double[] pars) {
        var newMolDataValues = new HashMap<String, MolDataValues>();

        for (var entry : molDataValues.entrySet()) {
            MolDataValues molData = entry.getValue();
            MolDataValues newMolData = new MolDataValues(molData.atom, molData.vector);
            newMolData.setTestModel(model);
            newMolDataValues.put(entry.getKey(), newMolData);
            MFModel testModel = newMolData.getTestModel();
            double[] resPars;
            if (useGlobalTau) {
                int nResPars = testModel.getNPars();
                resPars = new double[nResPars + 1];
                resPars[0] = globalTau;
                System.arraycopy(pars, 0, resPars, 1, nResPars);
            } else {
                resPars = pars;
            }

            for (var value : molData.getData()) {
                if (value instanceof R1R2NOEDataValue dValue) {
                    randomize(random, testModel, newMolData, dValue,resPars);
                } else {
                    var dValue = (DeuteriumDataValue) value;
                    randomize(random, testModel, newMolData, dValue,resPars);
                }
            }
        }
        return newMolDataValues;
    }

    private void randomize(Random random, MFModel testModel, MolDataValues newMolData, R1R2NOEDataValue dValue, double[] resPars) {
        RelaxEquations relaxObj = dValue.relaxObj;
        double[] valJ = testModel.calc(relaxObj.wValues, resPars);
        double r1 = relaxObj.R1(valJ);
        double rEx = testModel.includesEx() ? resPars[resPars.length - 1] : 0.0;
        double r2 = relaxObj.R2(valJ, rEx);
        double noe = relaxObj.NOE(valJ);
        dValue.randomize(newMolData, r1, r2, noe, random, 1.0);
    }

    private void randomize(Random random, MFModel testModel, MolDataValues newMolData, DeuteriumDataValue dValue, double[] resPars) {
        RelaxEquations relaxObj = dValue.relaxObj;
        double[] J = testModel.calc(relaxObj.wValues, resPars);
        double r1 = relaxObj.R1_D(J);
        double rEx = testModel.includesEx() ? resPars[resPars.length - 1] : 0.0;
        double r2 = relaxObj.R2_D(J);
        double rAP = relaxObj.Rap_D(J);
        double rQ = relaxObj.RQ_D(J);
        dValue.randomize(newMolData, r1, r2, rQ, rAP, random, 1.0);
    }

    public double value(double[] pars, double[][] values) {
        var score = score(pars, false);
        return score.value(getLambdaS(), getLambdaTau());
    }

    public double valueMultiResidue(double[] pars, double[][] values) {
        double sumSq = 0.0;
        int n = 0;
        boolean parsOK = true;
        int parStart = 1;
        for (MolDataValues molData : molDataValues.values()) {
            MFModel testModel = molData.getTestModel();
            int nResPars = testModel.getNPars();
            double[] resPars = new double[nResPars + 1];
            resPars[0] = pars[0];
            System.arraycopy(pars, parStart, resPars, 1, nResPars);
            parStart += nResPars;

            for (RelaxDataValue value : molData.getData()) {
                var dValue = (R1R2NOEDataValue) value;
                RelaxEquations relaxObj = dValue.relaxObj;
                double[] valJ = testModel.calc(relaxObj.wValues, resPars);
                double r1 = relaxObj.R1(valJ);
                double r2 = relaxObj.R2(valJ, 0.0);
                double noe = relaxObj.NOE(valJ);
                double delta2 = dValue.score2(r1, r2, noe);
                sumSq += delta2;
                n += 3;
            }
            if (!testModel.checkParConstraints()) {
                parsOK = false;
            }
        }
        double rms = n == 0 ? 0.0 : Math.sqrt(sumSq / n);
        if (!parsOK) {
            rms += n * 10.0;
        }
        return rms;
    }

    public double valueDMat(double[] pars, double[][] values) {
        if (diffusionType == OBLATE || diffusionType == PROLATE) {
            Arrays.sort(pars, 0, 2);
        } else if (diffusionType == ANISOTROPIC) {
            Arrays.sort(pars, 0, 3);
        }
        double[][] valD = null;
        double[][] valVT = null;
        if (diffusionType != DiffusionType.ISOTROPIC) {
            valD = parsToD(pars, diffusionType);
            valVT = parsToVT(pars, diffusionType);
        }
        double sumSq = 0.0;
        int nDiffPars = diffusionType.getNDiffusionPars() + diffusionType.getNAnglePars();
        int n = 0;
        for (MolDataValues molData : molDataValues.values()) {
            MFModel model = molData.getTestModel();
            for (RelaxDataValue value : molData.getData()) {
                R1R2NOEDataValue dValue  = (R1R2NOEDataValue) value;
                RelaxEquations relaxObj = dValue.relaxObj;
                double[] v = molData.vector;
                int nModelPars = model.getNPars();
                double[] resPars = new double[nModelPars + nDiffPars];
                System.arraycopy(pars, 0, resPars, 0, nDiffPars);
                resPars[resPars.length - 1] = 1.0; //Model 0: S2 = 1.0, all others null.
                double[] valJ;
                if (model instanceof MFModelIso) {
                    resPars[0] = 1.0 / (6.0 * resPars[0]);
                    valJ = model.calc(relaxObj.wValues, resPars);
                } else {
                    valJ = getJDiffusion(resPars, relaxObj, (MFModelAniso) model, v, valD, valVT);
                }
                double rhoExp = dValue.calcExpRho(valJ);
                double rhoPred = dValue.calcPredRho(valJ);
                double delta = rhoPred - rhoExp;
                sumSq += delta * delta;
                n++;
            }
        }

        return n == 0 ? 0.0 : Math.sqrt(sumSq / n);

    }

    public void dumpValues(double[] pars) {
        if (diffusionType == OBLATE || diffusionType == PROLATE) {
            Arrays.sort(pars, 0, 2);
        } else if (diffusionType == ANISOTROPIC) {
            Arrays.sort(pars, 0, 3);
        }
        double[][] valD = parsToD(pars, diffusionType);
        double[][] valVT = parsToVT(pars, diffusionType);
        int modelNum = 1;
        int nDiffPars = diffusionType.getNDiffusionPars() + diffusionType.getNAnglePars();
        for (MolDataValues molData : molDataValues.values()) {
            for (RelaxDataValue value : molData.getData()) {
                R1R2NOEDataValue dValue  = (R1R2NOEDataValue) value;
                RelaxEquations relaxObj = dValue.relaxObj;
                double[] v = molData.vector;
                double[] resPars = new double[nParsPerModel[modelNum] + nDiffPars];
                System.arraycopy(pars, 0, resPars, 0, nDiffPars);
                resPars[resPars.length - 1] = 1.0; //Model 0: S2 = 1.0, all others null.
                double[] valJ = getJDiffusion(resPars, relaxObj, modelNum, v, valD, valVT);
                double rhoExp = dValue.calcExpRho(valJ);
                double rhoPred = dValue.calcPredRho(valJ);
                System.out.println(rhoExp + " " + rhoPred + " " + (rhoExp - rhoPred) + " " + molData.specifier);
            }
        }
    }

    public double[] getPars() {
        return bestPars;
    }

    public double getAIC() {
        return bestAIC;
    }

    public double getChiSq() {
        return bestChiSq;
    }

    public double[] getParErrs() {
        return parErrs;
    }

    public Fitter getFitter() {
        return bestFitter;
    }

    private boolean checkParConstraint(double[] pars, int modelNum) {
        boolean result;
        double scale = 1.0;
        result = switch (modelNum) {
            case 2, 5 -> pars[2] < scale * pars[0];
            case 6 -> pars[2] < scale * pars[0] && pars[4] < scale * pars[0];
            default -> true;
        };
        return result;

    }

    public static double[][] guesses(int model, double tau) {
        double[] guess1 = {tau, 0.9};
        double[] lower1 = {tau / 10.0, 0.0};
        double[] upper1 = {tau * 10.0, 1.0};
        double[] guess2 = {tau, 0.9, tau / 40.0};
        double[] lower2 = {tau / 10.0, 0.0, tau / 1000.0};
        double[] upper2 = {tau * 10.0, 1.0, tau / 10.0};
        double[] guess5 = {tau, 0.9, tau / 40.0, 0.9};
        double[] lower5 = {tau / 10.0, 0.0, tau / 1000.0, 0.0};
        double[] upper5 = {tau * 10.0, 1.0, tau / 10.0, 1.0};
        double[] guess6 = {tau, 0.9, tau / 40.0, 0.9, tau / 40.0};
        double[] lower6 = {tau / 10.0, 0.0, tau / 1000.0, 0.0, tau / 1000.0};
        double[] upper6 = {tau * 10.0, 1.0, tau / 10.0, 1.0, tau / 10.0};
        double[][] guesses = {null, guess1, guess2, null, null, guess5, guess6};
        double[][] lower = {null, lower1, lower2, null, null, lower5, lower6};
        double[][] upper = {null, upper1, upper2, null, null, upper5, upper6};
        double[][] result = new double[3][];
        result[0] = guesses[model];
        result[1] = lower[model];
        result[2] = upper[model];
        return result;

    }

    public Optional<PointValuePair> fitResidueToModel(double[] start, double[] lower, double[] upper) {
        Fitter fitter = Fitter.getArrayFitter(this::value);
        try {
            return Optional.of(fitter.fit(start, lower, upper, 10.0));
        } catch (Exception ex) {
            ex.printStackTrace();
            return Optional.empty();
        }
    }

    public PointValuePair fitMultiResidueToModel(double[] start, double[] lower, double[] upper) {
        Fitter fitter = Fitter.getArrayFitter(this::valueMultiResidue);
        try {
            return fitter.fit(start, lower, upper, 10.0);
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }

    public PointValuePair fitDiffusion(double[] guesses) {
        Fitter fitter = Fitter.getArrayFitter(this::valueDMat);
        double[] lower = new double[guesses.length];
        double[] upper = new double[guesses.length];
        int nDiffPars = diffusionType.getNDiffusionPars();
        int nAnglePars = diffusionType.getNAnglePars();
        for (int i = 0; i < nDiffPars; i++) {
            lower[i] = guesses[i] / 2;
            upper[i] = guesses[i] * 2;
        }
        for (int i = nDiffPars; i < nDiffPars + nAnglePars; i++) {
            lower[i] = guesses[i] - Math.PI / 4.0;
            upper[i] = guesses[i] + Math.PI / 4.0;
        }
        try {
            PointValuePair result = fitter.fit(guesses, lower, upper, 10.0);
            bestPars = result.getPoint();
            bestChiSq = result.getValue();
            return result;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }
}
