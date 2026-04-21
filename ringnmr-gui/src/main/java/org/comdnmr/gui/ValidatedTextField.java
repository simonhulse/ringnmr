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

    private static final int DEFAULT_WIDTH = 60;
    private static final String VALID_STYLE = "";
    private static final String INVALID_STYLE = "-fx-text-box-border: red; -fx-focus-color: red; -fx-faint-focus-color: transparent;";

    private ValidationStrategy<T> validator;

    public ValidatedTextField(ValidationStrategy<T> validator) {
        this(validator, "");
    }

    public ValidatedTextField(ValidationStrategy<T> validator, String value) {
        this.validator = validator;
        setPrefWidth(DEFAULT_WIDTH);
        setMaxWidth(DEFAULT_WIDTH);
        setText(value);
        validateInput();
        bindListener();
    }

    private void bindListener() {
        textProperty().addListener(
            new ChangeListener<String>() {
                @Override
                public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
                    validateInput();
                }
            }
        );
    }

    public void validateInput() {
        boolean valid = getValue().isPresent();
        if (valid) {
            setStyle(VALID_STYLE);
        } else {
            setStyle(INVALID_STYLE);
        }
    }

    Optional<T> getValue() {
        return validator.getValue(getText());
    }
}
