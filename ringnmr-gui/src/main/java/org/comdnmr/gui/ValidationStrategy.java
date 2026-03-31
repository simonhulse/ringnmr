package org.comdnmr.gui;

import java.util.Optional;

public interface ValidationStrategy<T> {
    Optional<T> getValue(String input);
}
