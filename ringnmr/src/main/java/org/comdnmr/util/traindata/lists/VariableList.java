package org.comdnmr.util.traindata.lists;

import java.util.*;

public class VariableList {
    public String name;
    public ArrayList<ArrayList<Double>> values;

    public VariableList(String name) {
        this.name = name;
        this.values = new ArrayList<ArrayList<Double>>();
    }

    public void add(ArrayList<Double> value) {
        this.values.add(value);
    }

    public HashMap<String, Object> asHashMap() {
        HashMap<String, Object> map = new HashMap<>();
        map.put("name", this.name);
        map.put("values", this.values);
        return map;
    }
}
