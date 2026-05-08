package org.comdnmr.gui;

import java.util.Optional;

public class PositiveDecimalValidationStrategy implements ValidationStrategy<Double> {
    @Override
    public Optional<Double> getValue(String input) {
        try {
            Double value = Double.parseDouble(input);
            if (value >= 0) {
                return Optional.of(value);
            }
        } catch (NumberFormatException e) {}
        return Optional.empty();
    }
}
