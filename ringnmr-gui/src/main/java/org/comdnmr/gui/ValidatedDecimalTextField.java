package org.comdnmr.gui;

public class ValidatedDecimalTextField extends ValidatedTextField<Double> {

    public ValidatedDecimalTextField() {
        super(new PositiveDecimalValidationStrategy());
    }
}
