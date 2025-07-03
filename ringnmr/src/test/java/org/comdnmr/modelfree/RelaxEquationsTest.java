package org.comdnmr.modelfree;

import org.junit.Test;

public class RelaxEquationsTest {

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