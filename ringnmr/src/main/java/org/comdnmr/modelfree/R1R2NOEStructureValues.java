package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class R1R2NOEStructureValues extends StructureValues {

    @Override
    public Map<String, Double> estimateTau() {
        Map<Long, List<RelaxDataValue>> fieldMap = new HashMap<>();
        for (var entry : map.entrySet()) {
            for (RelaxDataValue rlxValue : entry.getValue().getData()) {
                long sfMHz = Math.round(rlxValue.relaxObj.getSF() / 1.0e6);
                fieldMap.computeIfAbsent(sfMHz, k -> new ArrayList<>()).add(rlxValue);
            }
        }
        int max = 0;
        long maxMHz = 0;
        for (long sfMHz : fieldMap.keySet()) {
            int size = fieldMap.get(sfMHz).size();
            if (size > max) {
                maxMHz = sfMHz;
                max = size;
            }
        }
        return max > 0
                ? CorrelationTime.estimateTau(maxMHz, "N", fieldMap.get(maxMHz))
                : Collections.emptyMap();
    }
}
