package org.comdnmr.modelfree;


import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class DeuteriumMapping {
    static Random random = new Random();
    static double[][] elements = {
            {0.0, 1.0, 4.0},
            {3.0 / 2.0, 5.0 / 2.0, 1.0},
            {0.0, 3.0, 0.0},
            {3.0 / 2.0, 1.0 / 2.0, 1.0}
    };

    private DeuteriumMapping() {

    }

    public static double[][] independentMapping(List<Double> rValueList, List<Double> errValueList, List<Double> fields) {
        int nRows = rValueList.size();
        int nCols = 3;
        double field = fields.get(0);
        double[] rValues = new double[nRows];
        for (int i = 0; i < rValues.length; i++) {
            rValues[i] = rValueList.get(i);
        }
        double[][] matrix = new double[nRows][nCols];
        for (int i = 0; i < nRows; i++) {
            double errScale = 1.0 / errValueList.get(i);
            for (int j = 0; j < nCols; j++) {
                matrix[i][j] = elements[i][j] * errScale;
            }
            rValues[i] *= errScale;
        }

        double scale = 3.0 * RelaxEquations.QCC2;

        OLSMultipleLinearRegression olsMultipleLinearRegression = new OLSMultipleLinearRegression();
        olsMultipleLinearRegression.setNoIntercept(true);
        olsMultipleLinearRegression.newSampleData(rValues, matrix);
        double[] jValues = olsMultipleLinearRegression.estimateRegressionParameters();
        double[] errors = olsMultipleLinearRegression.estimateRegressionParametersStandardErrors();
        double[] dFields = {0., field, field * 2};
        for (int i = 0; i < jValues.length; i++) {
            jValues[i] /= scale;
            errors[i] /= scale;
        }
        return new double[][]{dFields, jValues, errors};
    }

    public static double[][] jointMapping(List<Double> rValueList, List<Double> errValueList, List<Double> fields, boolean[] typeUsage) {
        int nRows = rValueList.size();
        double[] rValues = rValueList.stream()
                .mapToDouble(Double::doubleValue)
                .toArray();
        int nTypes = 0;
        for (boolean use : typeUsage) {
            if (use) {
                nTypes++;
            }
        }
        int[] typeMap = new int[nTypes];
        int k = 0;
        for (int i = 0; i < typeUsage.length; i++) {
            if (typeUsage[i]) {
                typeMap[k++] = i;
            }
        }

        int nFreqs = nRows / nTypes;

        List<Double> fieldList = new ArrayList<>();
        fieldList.add(0.0);
        int[] singleColumns = new int[nFreqs];
        int[] doubleColumns = new int[nFreqs];
        int iField = 0;
        for (var field : fields) {
            boolean matchSingle = false;
            boolean matchDouble = false;
            int jField = 0;
            for (var testField : fieldList) {
                if (Math.abs((field - testField) / testField) < 0.01) {
                    matchSingle = true;
                    singleColumns[iField] = jField;
                }
                if (Math.abs((field * 2.0 - testField) / testField) < 0.01) {
                    matchDouble = true;
                    doubleColumns[iField] = jField;
                }
                jField++;
            }
            if (!matchSingle) {
                singleColumns[iField] = fieldList.size();
                fieldList.add(field);
            }
            if (!matchDouble) {
                doubleColumns[iField] = fieldList.size();
                fieldList.add(field * 2.0);
            }

            iField++;
        }
        int nCols = fieldList.size();
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(nRows, nCols);

        for (int iFreq = 0; iFreq < nFreqs; iFreq++) {
            for (int iType = 0; iType < nTypes; iType++) {
                int row = iFreq * nTypes + iType;
                matrix.setEntry(row, 0, elements[typeMap[iType]][0]);
            }
        }
        for (int iFreq = 0; iFreq < nFreqs; iFreq++) {
            int singleColumn = singleColumns[iFreq];
            int doubleColumn = doubleColumns[iFreq];
            for (int iType = 0; iType < nTypes; iType++) {
                int row = iFreq * nTypes + iType;
                matrix.setEntry(row, singleColumn, elements[typeMap[iType]][1]);
            }
            for (int iType = 0; iType < nTypes; iType++) {
                int row = iFreq * nTypes + iType;
                matrix.setEntry(row, doubleColumn, elements[typeMap[iType]][2]);
            }
        }

        double scale = 3.0 * RelaxEquations.QCC2;

        for (int row = 0; row < nRows; row++) {
            double errScale = 1.0 / errValueList.get(row);
            for (int col = 0; col < nCols; col++) {
                matrix.multiplyEntry(row, col, errScale);
            }
        }

        try {
            int nReplicates = 500;
            double[][] jValuesRep = new double[nCols][nReplicates];
            double[][] errs = new double[nCols][nReplicates];
            double[] jValues = new double[nCols];
            OLSMultipleLinearRegression olsMultipleLinearRegression = new OLSMultipleLinearRegression();
            olsMultipleLinearRegression.setNoIntercept(true);
            for (int iRep = 0; iRep < nReplicates; iRep++) {
                double[] rTemp = new double[rValues.length];
                if (iRep == 0) {
                    for (int i = 0; i < nRows; i++) {
                        rTemp[i] = rValues[i] / errValueList.get(i);
                    }
                } else {
                    for (int row = 0; row < nRows; row++) {
                        rTemp[row] = (rValues[row] + random.nextGaussian() * errValueList.get(row)) / errValueList.get(row);

                    }
                }
                olsMultipleLinearRegression.newSampleData(rTemp, matrix.getData());
                double[] jValues1 = olsMultipleLinearRegression.estimateRegressionParameters();
                double[] errs1 = olsMultipleLinearRegression.estimateRegressionParametersStandardErrors();
                for (int i = 0; i < jValues1.length; i++) {
                    jValuesRep[i][iRep] = jValues1[i];
                    errs[i][iRep] = errs1[i];
                }
                if (iRep == 0) {
                    rValues = matrix.operate(jValues1);
                    for (int ii = 0; ii < rValues.length; ii++) {
                        rValues[ii] *= errValueList.get(ii);
                    }
                }
            }
            var jErrors = new double[jValues.length];
            var fitFields = new double[jValues.length];

            for (int i = 0; i < jValues.length; i++) {
                if (nReplicates == 1) {
                    jValues[i] = jValuesRep[i][0] / scale;
                    jErrors[i] = errs[i][0] / scale;
                } else {
                    DescriptiveStatistics sumStat = new DescriptiveStatistics(jValuesRep[i]);
                    jValues[i] = sumStat.getMean() / scale;
                    jErrors[i] = sumStat.getStandardDeviation() / scale;
                }
                fitFields[i] = fieldList.get(i);
            }
            double[][] jValuesOrig = new double[3][nFreqs * 3];
            for (int iFreq = 0; iFreq < nFreqs; iFreq++) {
                int singleColumn = singleColumns[iFreq];
                int doubleColumn = doubleColumns[iFreq];
                jValuesOrig[0][iFreq * 3] = fitFields[0];
                jValuesOrig[0][iFreq * 3 + 1] = fitFields[singleColumn];
                jValuesOrig[0][iFreq * 3 + 2] = fitFields[doubleColumn];
                jValuesOrig[1][iFreq * 3] = jValues[0];
                jValuesOrig[1][iFreq * 3 + 1] = jValues[singleColumn];
                jValuesOrig[1][iFreq * 3 + 2] = jValues[doubleColumn];
                jValuesOrig[2][iFreq * 3] = jErrors[0];
                jValuesOrig[2][iFreq * 3 + 1] = jErrors[singleColumn];
                jValuesOrig[2][iFreq * 3 + 2] = jErrors[doubleColumn];
            }
            double[] weights = new double[fitFields.length];
            Arrays.fill(weights, 1.0);
            return new double[][]{fitFields, jValues, jErrors, jValuesOrig[0], jValuesOrig[1], jValuesOrig[2], weights};
        } catch (IllegalArgumentException iAE) {
            return null;
        }

    }
}

