package org.comdnmr.gui;

import java.util.Optional;

public class PositiveIntegerValidationStrategy implements ValidationStrategy<Integer> {
    @Override
    public Optional<Integer> getValue(String input) {
        try {
            Integer value = Integer.parseInt(input);
            if (value > 0) {
                return Optional.of(value);
            }
        } catch (NumberFormatException e) {}
        return Optional.empty();
    }
}
