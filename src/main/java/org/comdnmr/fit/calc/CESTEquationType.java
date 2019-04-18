/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

/**
 *
 * @author Bruce Johnson
 */
public interface CESTEquationType extends EquationType {

    @Override
    public default double calculate(double[] par, int[] map, double[] X, int idNum, double field) {
        double[][] x = new double[3][1];
        //System.out.println(x.length + " " + x[0].length + " X " + X.length);
        x[0][0] = X[0];
        x[1][0] = X[1];
        x[2][0] = X[2];
        double[] y = calculate(par, map, x, idNum, field);
        return y[0];
    }

    @Override
    public default double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
        int nPars = CalcCEST.getNPars(map);
        double[] guesses = new double[nPars];
        for (int id = 0; id < map.length; id++) {
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            double[][] peaks = CESTEquations.cestPeakGuess(xy[0], xy[1], field);
            //                for (int i=0; i<peaks.length; i++) {
//                    for (int j=0; j<peaks[i].length; j++) {
//                        System.out.println("peaks " + i + " " + j + " = " + peaks[i][j]);
//                    }
//                }
            double tex = xValues[2][0];
            double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
            double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
            guesses[map1[0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
            guesses[map1[1]] = CESTEquations.cestPbGuess(peaks, yValues); //0.1; //pb
            guesses[map1[2]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
            guesses[map1[3]] = peaks[peaks.length - 1][0]; //400 * 2.0 * Math.PI; //deltaA
            guesses[map1[4]] = r1[0]; //2.4; //R1A
            guesses[map1[5]] = r1[1]; //2.4; //R1B
            guesses[map1[6]] = r2[0][0]; //20.0; //R2A
            guesses[map1[7]] = r2[1][0]; //100.0; //R2B
        }
//            for (int i=0; i<guesses.length; i++) {
//                System.out.println(guesses[i]);
//            }

        return guesses;
    }

    @Override
    public default double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
        double[][] boundaries = new double[2][guesses.length];
        for (int id = 0; id < map.length; id++) {
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            double[][] peaks = CESTEquations.cestPeakGuess(xy[0], xy[1], field);
            double dAbound = 0;
            double dBbound = 0;
            if (peaks.length > 1) {
                dAbound = (peaks[0][2] / field) / 2;
                dBbound = (peaks[1][2] / field) / 2;
            } else if (peaks.length == 1) {
                dAbound = (peaks[0][2] / field) / 2;
                dBbound = dAbound;
            }
            double tex = xValues[2][0];
            double r1A = guesses[map1[4]];
            double[] r1BouA = CESTEquations.r1Boundaries(r1A, tex, 0.1);
            double r1B = guesses[map1[5]];
            double[] r1BouB = CESTEquations.r1Boundaries(r1B, tex, 0.1);

            boundaries[0][map1[0]] = 1.0; //kex LB
            boundaries[1][map1[0]] = guesses[map1[0]] * 6; //kex UB
            boundaries[0][map1[1]] = 0.01; //pb LB
            boundaries[1][map1[1]] = 0.25; //pb UB //guesses[1] * 4;
            boundaries[0][map1[2]] = guesses[map1[2]] - dAbound; //deltaA LB
            boundaries[1][map1[2]] = guesses[map1[2]] + dAbound; //deltaA UB
            boundaries[0][map1[3]] = guesses[map1[3]] - dBbound; //deltaB LB
            boundaries[1][map1[3]] = guesses[map1[3]] + dBbound; //deltaB UB
            boundaries[0][map1[4]] = r1BouA[0]; //R1A LB
            boundaries[1][map1[4]] = r1BouA[1]; //R1A UB
            boundaries[0][map1[5]] = r1BouB[0]; //R1B LB
            boundaries[1][map1[5]] = r1BouB[1]; //R1B UB
            boundaries[0][map1[6]] = 2.0; //R2A LB
            boundaries[1][map1[6]] = guesses[map1[6]] * 6; //R2A UB
            boundaries[0][map1[7]] = 2.0; //R2B LB
            boundaries[1][map1[7]] = guesses[map1[7]] * 6; //R2B UB
            if (boundaries[1][map1[7]] < 200.0) {
                boundaries[1][map1[7]] = 200.0;
            }
        }

        return boundaries;
    }

    @Override
    public default double getRex(double[] pars, int[] map, double field) {
        return 0.0;
    }

    @Override
    public default double getKex(double[] pars) {
        return pars[0];
    }

    @Override
    public default double getKex(double[] pars, int id) {
        return pars[0];
    }

    @Override
    public default int[][] makeMap(int n) {
        int nP = 8;
        int[][] map = new int[n][nP];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < nP; j++) {
                map[i][j] = nP * i + j;
            }
        }
        return map;
    }

    @Override
    public default int[][] makeMap(int n, int m) {
        int nP = m;
        int[][] map = new int[n][nP];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < nP; j++) {
                map[i][0] = nP * i + j;
            }
        }
        return map;
    }

    @Override
    public default int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
        int[][] map = makeMap(1);
        return map;
    }
}
