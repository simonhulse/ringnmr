package org.comdnmr.gui;

import java.util.function.Predicate;

public class ValidatedDecimalTextField extends ValidatedTextField<Double> {

    public ValidatedDecimalTextField() {
        this("");
    }

    public ValidatedDecimalTextField(String initialValue) {
        super(new PositiveDecimalValidationStrategy(), initialValue);
    }

    public ValidatedDecimalTextField(Predicate<Double> additionalCheck, String initialValue) {
        super(new PositiveDecimalValidationStrategy(), additionalCheck, initialValue);
    }
}
