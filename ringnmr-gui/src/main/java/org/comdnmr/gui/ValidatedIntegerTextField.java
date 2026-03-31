package org.comdnmr.gui;

public class ValidatedIntegerTextField extends ValidatedTextField<Integer> {

    public ValidatedIntegerTextField() {
        super(new PositiveIntegerValidationStrategy());
    }
}
