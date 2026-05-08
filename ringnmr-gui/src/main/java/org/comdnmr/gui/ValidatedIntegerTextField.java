package org.comdnmr.gui;

import java.util.function.Predicate;

public class ValidatedIntegerTextField extends ValidatedTextField<Integer> {

    public ValidatedIntegerTextField() {
        this("");
    }

    public ValidatedIntegerTextField(String initialValue) {
        super(new PositiveIntegerValidationStrategy(), initialValue);
    }

    public ValidatedIntegerTextField(Predicate<Integer> additionalCheck, String initialValue) {
        super(new PositiveIntegerValidationStrategy(), additionalCheck, initialValue);
    }
}
