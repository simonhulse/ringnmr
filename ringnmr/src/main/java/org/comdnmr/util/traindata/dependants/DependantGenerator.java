package org.comdnmr.util.traindata.dependants;

public abstract class DependantGenerator {
    String name;

    final public String getName() { return name; }

    abstract public void updateState();

    abstract public double fetch(double value);
}
