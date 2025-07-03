package org.comdnmr.modelfree;

import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.Fitter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SpectralDensityCalculator {

    R1R2NOEDataValue relaxDataValue;

    public double value(double[] pars, double[][] values) {
        RelaxEquations relaxEquation = relaxDataValue.relaxObj;
        double ratio = RelaxEquations.GAMMA_N/ RelaxEquations.GAMMA_H;
        double f1 = 1.0 - ratio;
        double f2 = 1.0 + ratio;


        double jH = pars[2];
        double jHmN = jH  / Math.pow(f1,1.5);
        double jHpN = jH /  Math.pow(f2,1.5);

        double[] jValues = {pars[0], pars[1], jHmN, jH, jHpN};
        double r1 = relaxEquation.R1(jValues);
        double r2 = relaxEquation.R2(jValues, 0.0);
        double noe = relaxEquation.NOE(jValues);

        return relaxDataValue.score2(r1, r2, noe);
    }

    public PointValuePair fit(R1R2NOEDataValue dValue) {
        this.relaxDataValue = dValue;
        Fitter fitter = Fitter.getArrayFitter(this::value);
        double[][] jValues = calcJR1R2NOE(List.of(dValue));
        double j0 = jValues[1][0];
        double jN = jValues[1][2];
        double jH = jValues[1][1];
        double ratio = RelaxEquations.GAMMA_N/ RelaxEquations.GAMMA_H;
        double f1 = 1.0 - ratio;
        double f2 = 1.0 + ratio;

        double[] start = {j0, jN, jH};
        double[] lower = {j0 * 0.5, jN*0.5, jH*0.5};
        double[] upper = {j0 * 1.5, jN * 1.5, jH * 1.5};

        try {
            PointValuePair pointValuePair = fitter.fit(start, lower, upper, 10.0);
            double[] j = pointValuePair.getPoint();

            double jHr = j[2];
            double jHmN = jHr  / Math.pow(f1,1.5);
            double jHpN = jHr /  Math.pow(f2,1.5);

            double[] jResult = {j[0], j[1], jHmN, jHr, jHpN};
            return new PointValuePair(jResult, pointValuePair.getValue());
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }

    public static double[][] calcJR1R2NOE(List<RelaxDataValue> dataValues) {
        boolean useAverage = false;
        int nDataValues = dataValues.size();
        int stepSize = useAverage ? 2 : 3;
        int nFreq = useAverage ? 1 + stepSize * nDataValues : stepSize * nDataValues;
        double[][] result = new double[4][nFreq];
        int iField = 0;
        for (RelaxDataValue value : dataValues) {
            R1R2NOEDataValue relaxDataValue = (R1R2NOEDataValue) value;
            RelaxEquations relaxEq = relaxDataValue.relaxObj;
            double r1 = relaxDataValue.R1;
            double r2 = relaxDataValue.R2;
            double noe = relaxDataValue.NOE;
            double r1Err = relaxDataValue.R1err;
            double r2Err = relaxDataValue.R2err;
            double noeErr = relaxDataValue.NOEerr;
            double sigma = (noe - 1.0) * r1 * RelaxEquations.GAMMA_N / RelaxEquations.GAMMA_H;
            double sigmaErr = sigma * Math.sqrt(Math.pow((noeErr / (noe - 1.0)), 2) + Math.pow((r1Err / r1), 2));

            double d2 = relaxEq.getD2();
            double c2 = relaxEq.getC2();

            double j87H = 4.0 * sigma / (5.0 * d2);
            double j87Herr = 4.0 * sigmaErr / (5.0 * d2);

            double jNMul = 4.0 / (3.0 * d2 + 4.0 * c2);
            double jN = (r1 - 1.249 * sigma) * jNMul;
            double jNerr = jNMul * Math.sqrt(Math.pow(r1Err, 2) + Math.pow(1.249 * sigmaErr, 2));

            double j0Mul = 6.0 / (3.0 * d2 + 4.0 * c2);
            double j0 = j0Mul * (r2 - 0.5 * r1 - 0.454 * sigma);
            double j0Err = j0Mul * Math.sqrt(Math.pow(r2Err, 2) + Math.pow(0.5 * r1Err, 2) + Math.pow(0.454 * sigmaErr, 2));

            if (useAverage) {
                result[0][0] = 0.0;
                result[1][0] += j0;
                result[2][0] += j0Err * j0Err;
            } else {
                result[0][iField * stepSize] = 0.0;
                result[1][iField * stepSize] = j0;
                result[2][iField * stepSize] = j0Err;
            }

            result[0][iField * stepSize + 1] = 0.87 * relaxEq.getWI();
            result[1][iField * stepSize + 1] = j87H;
            result[2][iField * stepSize + 1] = j87Herr;

            result[0][iField * stepSize + 2] = relaxEq.getWS();
            result[1][iField * stepSize + 2] = jN;
            result[2][iField * stepSize + 2] = jNerr;
            iField++;
        }
        if (useAverage) {
            result[1][0] /= nDataValues;
            result[2][0] = Math.sqrt(result[2][0]);
        }
        Arrays.fill(result[3], 1.0);
        return result;
    }

    public static double[][] calcJDeuterium(List<RelaxDataValue> dataValues) {
        double[][] result = null;
        if (!dataValues.isEmpty()) {
            List<Double> rValues = new ArrayList<>();
            List<Double> errValues = new ArrayList<>();
            List<Double> fields = new ArrayList<>();
            boolean doIndependent = false;

            for (var value : dataValues) {
                var dValue = (DeuteriumDataValue) value;
                rValues.add(dValue.R1);
                rValues.add(dValue.R2);
                rValues.add(dValue.rQ);
                rValues.add(dValue.rAP);
                errValues.add(dValue.R1err);
                errValues.add(dValue.R2err);
                errValues.add(dValue.rQError);
                errValues.add(dValue.rAPError);
                fields.add(dValue.relaxObj.getSF() * RelaxEquations.GAMMA_D / RelaxEquations.GAMMA_H * 2.0 * Math.PI);
            }
            if (doIndependent && dataValues.size() == 1) {
                result = DeuteriumMapping.independentMapping(rValues, errValues, fields);
            } else {
                result = DeuteriumMapping.jointMapping(rValues, errValues, fields);
            }
        }
        return result;
    }
}
