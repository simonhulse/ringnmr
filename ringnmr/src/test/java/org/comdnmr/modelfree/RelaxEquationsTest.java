package org.comdnmr.modelfree;

import org.junit.Assert;
import org.junit.Test;

public class RelaxEquationsTest {
    double twoPI = 2.0 * Math.PI;
    double[] R1rhoCSATargets = {
        1.4148978688487408e-05, 7.221345551571136, 0.9586089899637809, 1.4148978688487293e-05,
        7.164909988657386, 1.8616125433308737, 1.414897868848692e-05, 6.985250038220472,
        193.56034510531762, 1.4148978688485414e-05, 6.292946892439036, 97.12185914367787
    };
    double[] R1rhoHNTargets = {
        0.0005196600618600594, 92.82942378966885, 12.322656458970256, 0.000519660061860058,
        92.1039606967329, 23.930510034543467, 0.0005196600618600531, 89.79448318660548,
        2488.163136700531, 0.0005196600618600338, 80.89512335997507, 1248.4738588553196
    };
    double[] R1rhoHHTargets = {
        0.00011720530140576303, 31.405449205080966, 4.168965929005294, 0.00011720530140576101,
        30.43014504038675, 455.87926068278733, 0.00011720530140575456, 27.462847707589887,
        229.54734968483226, 0.00011720530140572838, 18.700020481571983, 0.5590040134704212,
    };
    double[] R1rhoTaus = {1.0e-12, 1.0e-6, 1.0e-4};
    double[] R1rhoOmegaEs = {0.0, twoPI * 1.9e4, twoPI * 3.9e4, twoPI * 7.9e4};
    double R1rhoOmegaR = twoPI * 4.0e4;
    double R1rhoS2 = 0.9;
    double R1rhoSigmaH = 8.0 * 1.0e-6;  // Δσ
    double R1rhoWH = 6.0e8;

    @Test
    public void testR1rhoCSA() {
        RelaxEquations.setSigma("H", R1rhoSigmaH);
        RelaxEquations relax = new RelaxEquations(R1rhoWH, "H", "H");
        int i = 0;
        for (double R1rhoOmegaE: R1rhoOmegaEs) {
            for (double R1rhoTau: R1rhoTaus) {
                double target = R1rhoCSATargets[i];
                Assert.assertEquals(
                    target,
                    relax.r1RhoCSA(
                       R1rhoOmegaR, R1rhoOmegaE, R1rhoTau, R1rhoS2
                    ),
                    1.0e-6
                );
                i++;
            }
        }
    }

    @Test
    public void testR1rhoHN() {
        RelaxEquations relax = new RelaxEquations(R1rhoWH, "H", "N");
        int i = 0;
        for (double R1rhoOmegaE: R1rhoOmegaEs) {
            for (double R1rhoTau: R1rhoTaus) {
                double target = R1rhoHNTargets[i];
                Assert.assertEquals(
                    target,
                    relax.r1RhoIS(
                        R1rhoOmegaR, R1rhoOmegaE, R1rhoTau, R1rhoS2
                    ),
                    1.0e-6
                );
                i++;
            }
        }
    }

    @Test
    public void testR1rhoAA() {
        RelaxEquations.setSigma("H", R1rhoSigmaH);
        RelaxEquations relax = new RelaxEquations(R1rhoWH, "H", "H");
        int i = 0;
        for (double R1rhoOmegaE: R1rhoOmegaEs) {
            for (double R1rhoTau: R1rhoTaus) {
                double target = R1rhoHHTargets[i];
                Assert.assertEquals(
                        target,
                        relax.r1RhoAA(
                                R1rhoOmegaR, R1rhoOmegaE, R1rhoTau, R1rhoS2
                        ),
                        1.0e-6
                );
                i++;
            }
        }
    }

    @Test
    public void testNOE1() {
        RelaxEquations relaxEquations = new RelaxEquations(500.0e6, "H", "N");
        double noe = relaxEquations.NOE(5.0e-9);
        System.out.println(noe);
    }

    @Test
    public void testNOE2() {
        RelaxEquations relaxEquations = new RelaxEquations(500.0e6, "H", "N");
        double[] j = relaxEquations.getJ(5.0e-9);
        double noe = relaxEquations.NOE(j);
        System.out.println(noe);
    }
}
