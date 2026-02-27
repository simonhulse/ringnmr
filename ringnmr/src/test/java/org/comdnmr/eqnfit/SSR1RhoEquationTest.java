package org.comdnmr.eqnfit;

import java.util.Arrays;

import org.junit.Test;

public class SSR1RhoEquationTest {

    @Test
    public void testCalculate() {
        double[] taucs = new double[]{10.0e-9, 67.0e-6, 10.0e-3};  // 10ns, 67μs, 10ms
        double s2 = 0.3295;  // 109°
        double omegaR = 2.0 * Math.PI * 10.0e3;  // 10kHz
        double start = 3.0;
        double stop = 40.0;
        double step = 0.25;
        int n = (int) ((stop - start + 1) / step);
        double[] omega1s = new double[n];
        double[] omegaRs = new double[n];
        for (int i = 0; i < n; i++) {
            omega1s[i] = 2.0 * Math.PI * (start + i * step) * 1.0e3;
            omegaRs[i] = omegaR;
        }
        double[][] X = new double[4][n];
        X[0] = omega1s;
        X[3] = omegaRs;

        int[] map = new int[]{0, 1};
        System.out.printf("X[0]:%n%s%n", Arrays.toString(X[0]));
        for (double tauc : taucs) {
            double[] par = new double[]{tauc, s2};
            double[] r1RhoCSA = SSR1RhoEquation.CSA.calculate(par, map, X, 0);
            //System.out.printf("tauc: %s%n", tauc);
            //System.out.printf("%s%n%n", Arrays.toString(r1RhoCSA));
        }
    }
}
