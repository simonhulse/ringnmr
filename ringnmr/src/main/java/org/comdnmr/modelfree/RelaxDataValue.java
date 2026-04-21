/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

/**
 *
 * @author brucejohnson
 */
public abstract class RelaxDataValue  {

    final MolDataValues<?> molDataValue;
    double R1;
    double R1err;
    double R2;
    double R2err;
    final RelaxEquations relaxObj;

    public double getR1() { return R1; }

    public double getR1err() { return R1err; }

    public double getR2() { return R2; }

    public double getR2err() { return R2err; }

    public double getB0() { return relaxObj.getSF(); }

    public RelaxEquations getRelaxEquations() { return relaxObj; }

    public RelaxDataValue(MolDataValues<?> molDataValue, double r1,
                          double r1Error, double r2, double r2Error,
                          RelaxEquations relaxObj) {
        this.molDataValue = molDataValue;
        this.R1 = r1;
        this.R1err = r1Error;
        this.R2 = r2;
        this.R2err = r2Error;
        this.relaxObj = relaxObj;
    }


    public MolDataValues<?> getMolData() {
        return molDataValue;
    }

    public abstract double[] getObservables();
    public abstract double[] getObservableErrors();
    public abstract void setObservables(double[] values);

}
