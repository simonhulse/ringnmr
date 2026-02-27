package org.comdnmr.modelfree.models;

public class MFModelIso2sfx extends MFModelIso2sf {

    @Override
    public double[] calc(double[] omegas) {
        double[] J = new double[omegas.length];
        int j = 0;
        double ss2 = this.ss2;
        double sf2 = this.sf2 / sN;
        double s2 = ss2 * sf2;
        for (double omega : omegas) {
            omega *= 1.0e-9;
            double omega2 = omega * omega;
            double vM = s2 * tauM / (1.0 + omega2 * tauM * tauM);
            double vMS = sf2 * (1.0 - ss2) * (tauM * tauS * (tauM + tauS))
                    / (tauM * tauM * tauS * tauS * omega2 + (tauM + tauS) * (tauM + tauS));
            double vMF = (1.0 - sf2) * ss2 * (tauM * tauF * (tauM + tauF))
                    / (tauM * tauM * tauF * tauF * omega2 + (tauM + tauF) * (tauM + tauF));
            double tauMFS = tauF * (tauM + tauS) + tauM * tauS;
            double vMFS = (1.0 - sf2) * (1.0 - ss2) * (tauF * tauM * tauS * tauMFS)
                    / (tauF * tauF * tauM * tauM * tauS * tauS * omega2
                    + tauMFS * tauMFS);
            J[j++] = 0.4e-9 * (vM + vMS + vMF + vMFS);
        }
        complexityS = Math.abs(1.0 - sf2) + Math.abs(1.0 - ss2);
        complexityTau =
                Math.log10((tauS + tauPrime) / tauPrime) +
                        Math.log10((tauF + tauPrime) / tauPrime);
        return J;

    }
}
