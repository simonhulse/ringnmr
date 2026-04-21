package org.comdnmr.gui;

public class ValidatedIntegerTextField extends ValidatedTextField<Integer> {


    public ValidatedIntegerTextField() {
        this("");
    }

    public ValidatedIntegerTextField(String initialValue) {
        super(new PositiveIntegerValidationStrategy(), initialValue);
    }
}
