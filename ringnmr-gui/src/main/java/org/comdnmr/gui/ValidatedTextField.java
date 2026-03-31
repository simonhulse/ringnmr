package org.comdnmr.gui;

import java.util.Optional;

import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;

import javafx.scene.control.TextField;


@FunctionalInterface
interface InputValidator<T> {
    Optional<T> getValue(String input);
}

abstract class ValidatedTextField<T> extends TextField {

    private static final String INVALID_BORDER_COLOR = "red";

    private ValidationStrategy<T> validator;

    public ValidatedTextField(ValidationStrategy<T> validator) {
        this(validator, "");
    }

    public ValidatedTextField(ValidationStrategy<T> validator, String value) {
        this.validator = validator;

        textProperty().addListener(
            new ChangeListener<String>() {
                @Override
                public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
                    validateInput();
                }
            }
        );
        setText(value);
    }

    public void validateInput() {
        boolean valid = getValue().isPresent();
        if (valid || isDisabled()) {
            setStyle("");
        } else {
            setStyle(String.format("-fx-border-color: %s;", INVALID_BORDER_COLOR));
        }
    }

    Optional<T> getValue() {
        return validator.getValue(getText());
    }
}
