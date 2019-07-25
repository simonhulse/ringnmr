package org.comdnmr.fit.calc;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author brucejohnson
 */
public class RelaxEquations {

    public final static double MU0 = 4.0e-7 * Math.PI;
    public final static double GAMMA_N = -2.71e7;
    public final static double GAMMA_H = 2.68e8;
    public final static double PLANCK = 1.054e-34;
    public final static double R_HN = 1.01e-10;
    public final static double SIGMA = 165.0e-6;
    public final static Map<String, Double> GAMMA_MAP = new HashMap<>();
    public final static Map<String, Double> R_MAP = new HashMap<>();

    static {
        GAMMA_MAP.put("H", GAMMA_H);
        GAMMA_MAP.put("N", GAMMA_N);
        R_MAP.put("HN", R_HN);
        R_MAP.put("NH", R_HN);
    }

    private final double r;
    private final double d;
    private final double d2;
    private final double c;
    private final double c2;
    private final double gammaS;
    private final double gammaI;
    private final double sf;
    private final double wI;
    private final double wS;

    // don't use this yet.  Various inconsistencies with various different presentations of equations
    //   consider using scaled versions (smaller exponents)
    // add csa
    public RelaxEquations(double sf, String elem1, String elem2) {
        gammaI = GAMMA_MAP.get(elem1);
        gammaS = GAMMA_MAP.get(elem2);
        wI = sf * 2.0 * Math.PI;
        wS =wI * gammaS / gammaI;
        r = R_MAP.get(elem1 + elem2);
        d = MU0 * (gammaI * gammaS * PLANCK) / (4.0 * Math.PI * r * r * r);
        d2 = d * d;
        c = wS * SIGMA / Math.sqrt(3.0);
        c2 = c * c;

        this.sf = sf;
    }

    public double JModelFree(double w, double tau, double taui, double s2) {
        double value1 = s2 / (1.0 + w * w * taui * taui);
        double value2 = ((1.0 - s2) * (tau + taui) * tau) / ((tau + taui) * (tau + taui) + w * w * taui * taui * tau * tau);
        double value = 0.4 * taui * (value1 + value2);
        return value;
    }

    public double J(double w, double tau) {
        double value = 0.4 * tau / (1.0 + w * w * tau * tau);
        return value;
    }

    public double r2r1Ratio(double tau) {
        double num = 4.0 * J(0, tau) + J(wS - wI, tau) + 3.0 * J(wS, tau)
                + 6.0 * J(wI, tau) + 6.0 * J(wS + wI, tau)
                + (c2 / (3.0 * d2) * (4.0 * J(0, tau) + 3.0 * J(wS, tau)));

        double denom = 2.0 * J(wS - wI, tau) + 6.0 * J(wS, tau) + 12.0 * J(wS + wI, tau) + 2.0 * (c2 / (3.0 * d2) * J(wS, tau));
        return num / denom;
    }

    public double R1(double tau) {
        double dipolarContrib = d2 / 4.0 * (J(wI - wS, tau) + 3.0 * J(wS, tau) + 6.0 * J(wI + wS, tau));
        double csaContrib = c2 * J(wS, tau);
        return dipolarContrib + csaContrib;

    }

    public double R2(double tau) {
        double dipolarContrib = d2 / 8.0 * (4.0 * J(0.0, tau) + J(wI - wS, tau) + 3.0 * J(wS, tau)
                + 6.0 * J(wI, tau) + 6.0 * J(wI + wS, tau));
        double csaContrib = c2 / 6 * (4.0 * J(0.0, tau) + 3.0 * J(wS, tau));
        return dipolarContrib + csaContrib;
    }

}
