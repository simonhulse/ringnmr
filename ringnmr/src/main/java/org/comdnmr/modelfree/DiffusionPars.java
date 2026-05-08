/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import org.comdnmr.modelfree.RelaxFit.DiffusionType;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

/**
 *
 * @author brucejohnson
 */
public class DiffusionPars {

    DiffusionType diffType;
    public double[] dDiff;
    public double[] a;
    public double[] v;

    public DiffusionPars(DiffusionType diffType, double[] v) {
        this.diffType = diffType;
        this.v = v;
    }

    public DiffusionPars(DiffusionType diffType, double[][] D, double[][] VT, double[] v) {
        this.diffType = diffType;
        this.v = v;
        calcDiffusiond(diffType, D);
        calcDiffusiona(diffType, D, VT, v);
    }

    public void update(double[][] D, double[][] VT) {
        calcDiffusiond(diffType, D);
        calcDiffusiona(diffType, D, VT, v);

    }

    /**
     * Calculate the d array for the diffusion J(w) calculations. Equations from
     * the SI of Berlin K.; Longhini, A.; Dayie, T. K. and Fushman, D., J.
     * Biomol NMR, 2013.
     *
     * @param diffType DiffusionType. The type of diffusion: anisotropic,
     * prolate, or oblate.
     * @param D double[][]. The diffusion matrix.
     * @return double[]. The d array.
     */
    private void calcDiffusiond(DiffusionType diffType, double[][] D) {
        double Dx = D[0][0];
        double Dy = D[1][1];
        double Dz = D[2][2];
//        System.out.println("Dx = " + Dx + " Dy = " + Dy + " Dz = " + Dz);
        dDiff = new double[3];
        switch (diffType) {
            case ANISOTROPIC:
                dDiff = new double[5];
                double[] k = {Dy - Dx, Dz - Dx, (Dx + Dy + Dz) / 3.0, 0};
                k[3] = Math.sqrt(k[0] * k[0] - k[0] * k[1] + k[1] * k[1]);
                dDiff[0] = 4.0 * Dx + Dy + Dz;
                dDiff[1] = Dx + 4.0 * Dy + Dz;
                dDiff[2] = Dx + Dy + 4.0 * Dz;
                dDiff[3] = 6.0 * k[2] + 2.0 * k[3];
                dDiff[4] = 6.0 * k[2] - 2.0 * k[3];
                break;
            case PROLATE: {
                //Dxx = Dyy
                double Dpar = Dz;
                double Dperp = Dx;
                dDiff[0] = 5.0 * Dperp + Dpar;
                dDiff[1] = 2.0 * Dperp + 4.0 * Dpar;
                dDiff[2] = 6.0 * Dperp;
                break;
            }
            case OBLATE: {
                //Dyy = Dzz
                double Dpar = Dx;
                double Dperp = Dz;
                dDiff[0] = 5.0 * Dperp + Dpar;
                dDiff[1] = 2.0 * Dperp + 4.0 * Dpar;
                dDiff[2] = 6.0 * Dperp;
                break;
            }
            default:
                break;
        }
    }

    /**
     * Calculate the a array for the diffusion J(w) calculations. Equations from
     * the SI of Berlin K.; Longhini, A.; Dayie, T. K. and Fushman, D., J.
     * Biomol NMR, 2013.
     *
     * @param diffType DiffusionType. The type of diffusion: anisotropic,
     * prolate, or oblate.
     * @param D double[][]. The diffusion matrix.
     * @param vec double[]. The unit vector for the SI bond.
     * @return double[]. The a array.
     */
    private void calcDiffusiona(DiffusionType diffType, double[][] D, double[][] VT, double[] vec) {
        double Dx = D[0][0];
        double Dy = D[1][1];
        double Dz = D[2][2];
        DMatrixRMaj vec1 = new DMatrixRMaj(vec);
        DMatrixRMaj VT1 = new DMatrixRMaj(VT);
        DMatrixRMaj v = new DMatrixRMaj(vec1.numRows, vec1.numCols);
        CommonOps_DDRM.mult(VT1, vec1, v);
        double vx = v.get(0, 0); //vec[0]; //vec.getX();
        double vy = v.get(1, 0); //vec[1]; //vec.getY();
        double vz = v.get(2, 0); //vec[2]; //vec.getZ();
        double vx2 = vx * vx;
        double vy2 = vy * vy;
        double vz2 = vz * vz;
        a = new double[3];
        switch (diffType) {
            case ANISOTROPIC:
                a = new double[5];
                double[] k = {Dy - Dx, Dz - Dx, (Dx + Dy + Dz) / 3.0, 0.0};
                k[3] = Math.sqrt(k[0] * k[0] - k[0] * k[1] + k[1] * k[1]);
                double[] delta = {(-k[0] - k[1]) / k[3], (2.0 * k[0] - k[1]) / k[3], (2.0 * k[1] - k[0]) / k[3]};
                if (k[0] <= 1e-12 && k[1] <= 1e-12 && k[3] <= 1e-12) {
                    delta = new double[3];
                }
//                for (int i=0; i<delta.length; i++) {
//                    System.out.print("del " + i + " = " + delta[i] + " ");
//                }
//                System.out.println();
                a[0] = 3.0 * (vy2) * (vz2);
                a[1] = 3.0 * (vx2) * (vz2);
                a[2] = 3.0 * (vx2) * (vy2);
                double p1 = 0.25 * (3.0 * ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2)) - 1.0);
                double val1 = delta[0] * (3 * vx2 * vx2 + 2 * a[0] - 1.0);
                double val2 = delta[1] * (3 * vy2 * vy2 + 2 * a[1] - 1.0);
                double val3 = delta[2] * (3 * vz2 * vz2 + 2 * a[2] - 1.0);
                double p2 = (1.0 / 12.0) * (val1 + val2 + val3);
                a[3] = p1 - p2;
                a[4] = p1 + p2;
                break;
            case PROLATE:
                //Dxx = Dyy
                a[0] = 3.0 * (vz2) * (1.0 - vz2);
                a[1] = 0.75 * (1.0 - vz2) * (1.0 - vz2);
                a[2] = 0.25 * (3.0 * vz2 - 1.0) * (3.0 * vz2 - 1.0);
                break;
            case OBLATE:
                //Dyy = Dzz
                a[0] = 3.0 * (vx2) * (1 - vx2);
                a[1] = 0.75 * (1 - vx2) * (1 - vx2);
                a[2] = 0.25 * (3 * vx2 - 1) * (3 * vx2 - 1);
                break;
            default:
                break;
        }
//        for (int i=0; i<a.length; i++) {
//            System.out.println("a " + i + ": " + a[i]);
//        }
    }

    /**
     * Calculate the e array for the diffusion J(w) calculations. Equations from
     * the SI of Berlin K.; Longhini, A.; Dayie, T. K. and Fushman, D., J.
     * Biomol NMR, 2013.
     *
     * @param tauLoc double. The internal correlation time.
     * @return double[]. The e array.
     */
    public double[] calcDiffusione(double tauLoc) {
        double[] e = new double[dDiff.length];
        for (int i = 0; i < e.length; i++) {
            e[i] = tauLoc / (dDiff[i] * tauLoc + 1.0);
        }
        return e;
    }

    public double[] getDf(double w2) {
        double[] Df = new double[dDiff.length];
        for (int d = 0; d < Df.length; d++) {
            if (w2 > 0.0) {
                Df[d] = dDiff[d] / (dDiff[d] * dDiff[d] + w2);
            } else {
                Df[d] = 1.0 / dDiff[d];
            }
        }
        return Df;
    }

}
