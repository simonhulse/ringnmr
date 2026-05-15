package org.comdnmr.modelfree;

import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.Fitter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SpectralDensityCalculator {
    static boolean useR1 = true;
    static boolean useR2 = true;
    static boolean useNOE = true;
    static boolean useRQ = true;
    static boolean useRAP = true;

    R1R2NOEDataValue relaxDataValue;

    public static void setUseRQ(boolean state) {
        useRQ = state;
    }

    public double value(double[] pars, double[][] values) {
        RelaxEquations relaxEquation = relaxDataValue.relaxObj;
        double ratio = RelaxEquations.GAMMA_N / RelaxEquations.GAMMA_H;
        double f1 = 1.0 - ratio;
        double f2 = 1.0 + ratio;


        double jH = pars[2];
        double jHmN = jH / Math.pow(f1, 1.5);
        double jHpN = jH / Math.pow(f2, 1.5);

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
        double ratio = RelaxEquations.GAMMA_N / RelaxEquations.GAMMA_H;
        double f1 = 1.0 - ratio;
        double f2 = 1.0 + ratio;

        double[] start = {j0, jN, jH};
        double[] lower = {j0 * 0.5, jN * 0.5, jH * 0.5};
        double[] upper = {j0 * 1.5, jN * 1.5, jH * 1.5};

        try {
            PointValuePair pointValuePair = fitter.fit(start, lower, upper, 10.0);
            double[] j = pointValuePair.getPoint();

            double jHr = j[2];
            double jHmN = jHr / Math.pow(f1, 1.5);
            double jHpN = jHr / Math.pow(f2, 1.5);

            double[] jResult = {j[0], j[1], jHmN, jHr, jHpN};
            return new PointValuePair(jResult, pointValuePair.getValue());
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }

    public static double[][] calcJR1R2NOE(List<R1R2NOEDataValue> dataValues) {
        int nFields = dataValues.size();
        double[][] result = new double[4][3 * nFields];
        for (int iField = 0; iField < nFields; iField++) {
            R1R2NOEDataValue relaxDataValue = dataValues.get(iField);
            RelaxEquations relaxEq = relaxDataValue.relaxObj;
            double r1 = relaxDataValue.R1;
            double r2 = relaxDataValue.R2;
            double noe = relaxDataValue.NOE;
            double r1Err = relaxDataValue.R1err;
            double r2Err = relaxDataValue.R2err;
            double noeErr = relaxDataValue.NOEerr;

            double sigma    = (noe - 1.0) * r1 * RelaxEquations.GAMMA_N / RelaxEquations.GAMMA_H;
            double sigmaErr = sigma * Math.sqrt(sq(noeErr / (noe - 1.0)) + sq(r1Err / r1));

            double d2 = relaxEq.getD2();
            double c2 = relaxEq.getC2();

            double j0Mul  = 6.0 / (3.0 * d2 + 4.0 * c2);
            double j0     = j0Mul * (r2 - 0.5 * r1 - 0.454 * sigma);
            double j0Err  = j0Mul * Math.sqrt(sq(r2Err) + sq(0.5 * r1Err) + sq(0.454 * sigmaErr));

            double j87H    = 4.0 * sigma / (5.0 * d2);
            double j87Herr = 4.0 * sigmaErr / (5.0 * d2);

            double jNMul = 4.0 / (3.0 * d2 + 4.0 * c2);
            double jN    = (r1 - 1.249 * sigma) * jNMul;
            double jNerr = jNMul * Math.sqrt(sq(r1Err) + sq(1.249 * sigmaErr));

            result[0][iField * 3]     = 0.0;
            result[1][iField * 3]     = j0;
            result[2][iField * 3]     = j0Err;

            result[0][iField * 3 + 1] = 0.87 * relaxEq.getWI();
            result[1][iField * 3 + 1] = j87H;
            result[2][iField * 3 + 1] = j87Herr;

            result[0][iField * 3 + 2] = relaxEq.getWS();
            result[1][iField * 3 + 2] = jN;
            result[2][iField * 3 + 2] = jNerr;
        }
        Arrays.fill(result[3], 1.0);
        return result;
    }

    public static int getNData(List<? extends RelaxDataValue> dataValues) {
        int n = 0;
        if (!dataValues.isEmpty()) {
            boolean[] typeUsage;
            if (dataValues.get(0) instanceof DeuteriumDataValue) {
                typeUsage = new boolean[]{useR1, useR2, useRQ, useRAP};
            } else {
                typeUsage = new boolean[]{useR1, useR2, useNOE};
            }
            int nActive = 0;
            for (boolean type : typeUsage) {
                if (type) {
                    nActive++;
                }
            }
            n = dataValues.size() * nActive;
        }
        return n;
    }

    /**
     * Computes J(ω) values using the averaged-J(0) strategy: a single J(0) is
     * derived from a no-intercept weighted least-squares regression of the
     * Γ = R2 − 0.5·R1 − 0.454·σ values across all fields, yielding a
     * J-value vector of length 1 + 2F (F = number of fields).
     *
     * <p>Layout: {@code [J(0), J(0.87ωH₀), J(ωN₀), J(0.87ωH₁), J(ωN₁), ...]}.
     *
     * <p>J(0) uncertainty is estimated by jackknife (leave-one-field-out) using
     * precision weights 1/σ²(Γᵢ) regardless of {@code j0Weights}, so the
     * reported error reflects the spread across fields rather than the
     * particular bootstrap draw. For F = 1, error propagation is used instead.
     *
     * @param dataValues the per-field relaxation observables; must not be empty
     * @param j0Weights  per-field weights applied in the J(0) regression
     *                   (bootstrap selection counts or Dirichlet draws);
     *                   {@code null} means uniform weights
     * @return {@code double[4][1+2F]}: rows are frequencies, J values, errors,
     *         and initial weights (all 1.0; overwritten by chi-sq weights later)
     */
    public static double[][] calcJR1R2NOEAveraged(List<R1R2NOEDataValue> dataValues, double[] j0Weights) {
        int nFields = dataValues.size();
        double[][] result = new double[4][1 + 2 * nFields];

        double[] gamma    = new double[nFields];
        double[] gammaErr = new double[nFields];
        double[] s        = new double[nFields];

        for (int i = 0; i < nFields; i++) {
            R1R2NOEDataValue dv = dataValues.get(i);
            RelaxEquations relaxEq = dv.relaxObj;
            double r1 = dv.R1, r1Err = dv.R1err;
            double r2 = dv.R2, r2Err = dv.R2err;
            double noe = dv.NOE, noeErr = dv.NOEerr;

            double sigma    = (noe - 1.0) * r1 * RelaxEquations.GAMMA_N / RelaxEquations.GAMMA_H;
            double sigmaErr = sigma * Math.sqrt(
                sq(noeErr / (noe - 1.0)) + sq(r1Err / r1)
            );

            double d2     = relaxEq.getD2();
            double c2     = relaxEq.getC2();
            double j0Mul  = 6.0 / (3.0 * d2 + 4.0 * c2);

            gamma[i]    = r2 - 0.5 * r1 - 0.454 * sigma;
            gammaErr[i] = Math.sqrt(sq(r2Err) + sq(0.5 * r1Err) + sq(0.454 * sigmaErr));
            s[i]        = 1.0 / j0Mul;
        }

        // No-intercept WLS: J(0) = Σ(pᵢ·sᵢ·Γᵢ) / Σ(pᵢ·sᵢ²), pᵢ = w[i]/σ²(Γᵢ)
        double num = 0.0, den = 0.0;
        for (int i = 0; i < nFields; i++) {
            double p = (j0Weights != null ? j0Weights[i] : 1.0) / sq(gammaErr[i]);
            num += p * s[i] * gamma[i];
            den += p * sq(s[i]);
        }
        double j0 = den > 0.0 ? num / den : 0.0;

        // Jackknife error always uses precision-only weights for consistency
        double j0Err;
        if (nFields == 1) {
            j0Err = gammaErr[0] / s[0];
        } else {
            double sumSq = 0.0;
            for (int k = 0; k < nFields; k++) {
                double numK = 0.0, denK = 0.0;
                for (int i = 0; i < nFields; i++) {
                    if (i == k) continue;
                    double p = 1.0 / sq(gammaErr[i]);
                    numK += p * s[i] * gamma[i];
                    denK += p * sq(s[i]);
                }
                double j0k = numK / denK;
                sumSq += sq(j0k - j0);
            }
            j0Err = Math.sqrt((double) (nFields - 1) / nFields * sumSq);
        }

        result[0][0] = 0.0;
        result[1][0] = j0;
        result[2][0] = j0Err;
        result[3][0] = 1.0;

        for (int i = 0; i < nFields; i++) {
            R1R2NOEDataValue dv = dataValues.get(i);
            RelaxEquations relaxEq = dv.relaxObj;
            double r1 = dv.R1, r1Err = dv.R1err;
            double noe = dv.NOE, noeErr = dv.NOEerr;

            double sigma    = (noe - 1.0) * r1 * RelaxEquations.GAMMA_N / RelaxEquations.GAMMA_H;
            double sigmaErr = sigma * Math.sqrt(
                sq(noeErr / (noe - 1.0)) + sq(r1Err / r1)
            );

            double d2 = relaxEq.getD2();
            double c2 = relaxEq.getC2();

            double j87H    = 4.0 * sigma / (5.0 * d2);
            double j87Herr = 4.0 * Math.abs(sigmaErr) / (5.0 * d2);

            double jNMul = 4.0 / (3.0 * d2 + 4.0 * c2);
            double jN    = (r1 - 1.249 * sigma) * jNMul;
            double jNerr = jNMul * Math.sqrt(sq(r1Err) + sq(1.249 * sigmaErr));

            result[0][1 + i * 2]     = 0.87 * relaxEq.getWI();
            result[1][1 + i * 2]     = j87H;
            result[2][1 + i * 2]     = j87Herr;
            result[3][1 + i * 2]     = 1.0;

            result[0][1 + i * 2 + 1] = relaxEq.getWS();
            result[1][1 + i * 2 + 1] = jN;
            result[2][1 + i * 2 + 1] = jNerr;
            result[3][1 + i * 2 + 1] = 1.0;
        }

        return result;
    }

    private static double sq(double x) { return x * x; }

    public static double[][] calcJDeuterium(List<DeuteriumDataValue> dataValues) {
        double[][] result = null;
        if (!dataValues.isEmpty()) {
            List<Double> rValues = new ArrayList<>();
            List<Double> errValues = new ArrayList<>();
            List<Double> fields = new ArrayList<>();
            boolean doIndependent = false;
            boolean[] typeUsage = {useR1, useR2, useRQ, useRAP};

            for (var dValue : dataValues) {
                if (typeUsage[0]) {
                    rValues.add(dValue.R1);
                    errValues.add(dValue.R1err);
                }
                if (typeUsage[1]) {
                    rValues.add(dValue.R2);
                    errValues.add(dValue.R2err);
                }

                if (typeUsage[2]) {
                    rValues.add(dValue.rQ);
                    errValues.add(dValue.rQError);
                }
                if (typeUsage[3]) {
                    rValues.add(dValue.rAP);
                    errValues.add(dValue.rAPError);
                }
                fields.add(dValue.relaxObj.getSF() * RelaxEquations.GAMMA_D / RelaxEquations.GAMMA_H * 2.0 * Math.PI);
            }
            if (doIndependent || dataValues.size() == 1) {
                result = DeuteriumMapping.independentMapping(rValues, errValues, fields);
            } else {
                result = DeuteriumMapping.jointMapping(rValues, errValues, fields, typeUsage);
            }
        }
        return result;
    }
}
