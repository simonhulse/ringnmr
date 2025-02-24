package org.comdnmr.util.traindata;

import java.util.*;

public class DataList<T> {
    public String name;
    public ArrayList<T> values;

    public DataList(String name) {
        this.name = name;
        this.values = new ArrayList<T>();
    }

    public void add(T value) {
        this.values.add(value);
    }

    public HashMap<String, Object> asHashMap() {
        HashMap<String, Object> map = new HashMap<>();
        map.put("name", this.name);
        map.put("values", this.values);
        return map;
    }
}
