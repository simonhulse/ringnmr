package org.comdnmr.util.traindata.lists;

import java.util.*;

public class ProfileList {
    public final String name = "profile";
    public ArrayList<double[][]> values;

    public ProfileList() {
        this.values = new ArrayList<>();
    }

    public void add(double[][] value) {
        this.values.add(value);
    }

    public HashMap<String, Object> asHashMap() {
        HashMap<String, Object> map = new HashMap<>();
        map.put("name", this.name);
        map.put("values", this.values);
        return map;
    }
}
