package org.comdnmr.gui;

public class ValidatedDecimalTextField extends ValidatedTextField<Double> {

    public ValidatedDecimalTextField() {
        this("");
    }

    public ValidatedDecimalTextField(String initialValue) {
        super(new PositiveDecimalValidationStrategy(), initialValue);
    }
}
