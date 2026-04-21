package org.comdnmr.modelfree;

import java.util.Map;

public class DeuteriumStructureValues extends StructureValues {

    @Override
    public Map<String, Double> estimateTau() {
        return Map.of("tau", 10.0);
    }
}
