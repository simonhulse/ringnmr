package org.comdnmr.util.traindata.lists;

import java.util.*;

public class ParameterList {
    public String name;
    public ArrayList<Double> values;

    public ParameterList(String name) {
        this.name = name;
        this.values = new ArrayList<Double>();
    }

    public void add(double value) {
        this.values.add(value);
    }

    public HashMap<String, Object> asHashMap() {
        HashMap<String, Object> map = new HashMap<String, Object>();
        map.put("name", this.name);
        map.put("values", this.values);
        return map;
    }
}
