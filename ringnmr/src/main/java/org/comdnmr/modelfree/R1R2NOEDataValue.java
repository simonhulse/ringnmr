package org.comdnmr.modelfree;

import java.util.Random;

import static org.comdnmr.modelfree.RelaxEquations.*;

public class R1R2NOEDataValue extends RelaxDataValue {
    double NOE;
    final double NOEerr;
    double rhoError = 0.0;

    public R1R2NOEDataValue(MolDataValues molDataValue, double r1,
                          double r1Error, double r2, double r2Error, double noe, double noeError,
                          RelaxEquations relaxObj) {
        super(molDataValue, r1,r1Error,r2,r2Error, relaxObj);
        this.NOE = noe;
        this.NOEerr = noeError;
    }

    public double getNOE() { return NOE; }

    public double getNOEerr() { return NOEerr; }

    public void randomize(MolDataValues molData, double r1, double r2,
                          double noe, Random random, double scale) {
        double newR1 = r1 + random.nextGaussian() * scale * R1err;
        double newR2 = r2 + random.nextGaussian() * scale * R2err;
        double newNOE = noe + random.nextGaussian() * scale * NOEerr;
        RelaxDataValue newValue = new R1R2NOEDataValue(molData, newR1, R1err,
                newR2, R2err, newNOE, NOEerr, relaxObj);
        molData.addData(newValue);
    }

    @Override
    public String toString() {
        return String.format("R1 %.2f +/- %.2f R2 %.2f +/- %.2f NOE %.2f +/- %.2f",
                R1, R1err, R2, R2err, NOE, NOEerr);
    }

    public String deltaString(double r1P, double r2P, double noeP) {
        return String.format("R1 %.2f +/- %.2f R2 %.2f +/- %.2f NOE %.2f +/- %.2f",
                r1P, (r1P - R1), r2P, (r2P - R2), noeP, (noeP - NOE));
    }

    public void updateRhoError(double[] J) {
        double gammaS = relaxObj.getGammaS();
        double gammaI = relaxObj.getGammaI();
        double w1 = (J[ImS] + 6 * J[IpS]) / (6 * J[IpS] - J[ImS]);
        double w2 = (6 * J[I]) / (6 * J[IpS] - J[ImS]);
        double f = (gammaS / gammaI) * R1 * (NOE - 1);
        double rhoExp1 = 2 * R2 - R1 - w2 * f;
        double rhoExp2 = R1 - w1 * f;
        double rhoExpErr1 = Math.sqrt(2 * 2 * R2err * R2err + R1err * R1err
                + w2 * w2 * (f * f * ((R1err / R1) * (R1err / R1) + (NOEerr / NOE) * (NOEerr / NOE)) + R1err * R1err));
        double rhoExpErr2 = Math.sqrt(R1err * R1err
                + w1 * w1 * (f * f * ((R1err / R1) * (R1err / R1) + (NOEerr / NOE) * (NOEerr / NOE)) + R1err * R1err));
        double rhoExp = rhoExp1 / rhoExp2;
        rhoError = Math.sqrt(rhoExp * rhoExp * ((rhoExpErr1 / rhoExp1) * (rhoExpErr1 / rhoExp1) + (rhoExpErr2 / rhoExp2) * (rhoExpErr2 / rhoExp2)));

    }

    public double calcExpRho(double[] J) {
        return relaxObj.calcRhoExp(R1, R2, NOE, J);
    }

    public double calcPredRho(double[] J) {
        return relaxObj.calcRhoPred(J);
    }

    public double score2(double r1P, double r2P, double noeP) {
        double r1D = (r1P - R1) / R1err;
        double r2D = (r2P - R2) / R2err;
        double noeD = (noeP - NOE) / NOEerr;
        return r1D * r1D + r2D * r2D + noeD * noeD;
    }

    public double score2(double r1P, double r2P, double noeP, double[] counts) {
        double r1D = (r1P - R1) / R1err;
        double r2D = (r2P - R2) / R2err;
        double noeD = (noeP - NOE) / NOEerr;
        return counts[0] * r1D * r1D + counts[1] * r2D * r2D + counts[2] * noeD * noeD;
    }

    public double scoreAbs(double r1P, double r2P, double noeP) {
        double r1D = (r1P - R1) / R1err;
        double r2D = (r2P - R2) / R2err;
        double noeD = (noeP - NOE) / NOEerr;
        return Math.abs(r1D) + Math.abs(r2D) + Math.abs(noeD);
    }
}
