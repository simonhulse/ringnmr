package org.comdnmr.modelfree;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public abstract class StructureValues {

    protected final Map<String, MolDataValues<? extends RelaxDataValue>> map = new HashMap<>();

    public boolean isEmpty() { return map.isEmpty(); }
    public void clear() { map.clear(); }
    public boolean containsKey(String key) { return map.containsKey(key); }
    public MolDataValues<? extends RelaxDataValue> get(String key) { return map.get(key); }
    public void put(String key, MolDataValues<? extends RelaxDataValue> value) { map.put(key, value); }
    public int size() { return map.size(); }
    public Set<Map.Entry<String, MolDataValues<? extends RelaxDataValue>>> entrySet() { return map.entrySet(); }
    public Collection<MolDataValues<? extends RelaxDataValue>> values() { return map.values(); }

    public abstract Map<String, Double> estimateTau();
}
