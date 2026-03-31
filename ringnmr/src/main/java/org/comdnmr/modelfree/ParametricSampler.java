package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ParametricSampler extends BootstrapSampler {

    private int nGaussianCalls = 0;
    private double cachedZ1 = 0.0;

    private Class<? extends RelaxDataValue> dataType;
    private Map<String, List<Double>> originalDataMap;

    public ParametricSampler(MolDataValues data) {
        super(data);
        setOriginalData();
    }

    private void setOriginalData() {
        try {
            dataType = data.dataValues.get(0).getClass();
        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("`MolDataValues` object has no dataValues!");
        }

        originalDataMap = new HashMap<>();

        originalDataMap.put("R1", new ArrayList<>());
        originalDataMap.put("R2", new ArrayList<>());
        for (RelaxDataValue dataValue : data.dataValues) {
            originalDataMap.get("R1").add(dataValue.R1);
            originalDataMap.get("R2").add(dataValue.R2);
        }

        if (dataType.equals(R1R2NOEDataValue.class)) {
            originalDataMap.put("NOE", new ArrayList<>());
            for (RelaxDataValue dataValue : data.dataValues) {
                R1R2NOEDataValue value = (R1R2NOEDataValue) dataValue;
                originalDataMap.get("NOE").add(value.NOE);
            }
        } else if (dataType.equals(DeuteriumDataValue.class)) {
            originalDataMap.put("rQ", new ArrayList<>());
            originalDataMap.put("rAP", new ArrayList<>());
            for (RelaxDataValue dataValue : data.dataValues) {
                DeuteriumDataValue value = (DeuteriumDataValue) dataValue;
                originalDataMap.get("rAP").add(value.rAP);
                originalDataMap.get("rQ").add(value.rQ);
            }
        } else {
            throw new AssertionError("Unexpected subclass of `RelaxDataValue` detected!");
        }
    }


    // Compute two independed normally distributed random values using the
    // Box-Muller transform
    // z1 is cached, so computation of z0 and z1 is required with every
    // odd-numbered call
    private double nextGaussian(double mean, double stdev) {
        double z;
        if (nGaussianCalls % 2 == 1) {
            z = cachedZ1;
        } else {
            double u1 = rng.nextDouble();
            double u2 = rng.nextDouble();
            z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
            cachedZ1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);
        }
        nGaussianCalls++;
        return mean + stdev * z;
    }

    public MolDataValues sample() {
        List<RelaxDataValue> dataValues = data.getData();
        for (int i = 0; i < getNFields(); i++) {
            RelaxDataValue dataValue = dataValues.get(i);
            dataValue.R1 = nextGaussian(originalDataMap.get("R1").get(i), dataValue.R1err);
            dataValue.R2 = nextGaussian(originalDataMap.get("R2").get(i), dataValue.R2err);

            if (dataType.equals(R1R2NOEDataValue.class)) {
                R1R2NOEDataValue value = (R1R2NOEDataValue) dataValue;
                value.NOE = nextGaussian(originalDataMap.get("NOE").get(i), value.NOEerr);
            } else if (dataType.equals(DeuteriumDataValue.class)) {
                DeuteriumDataValue value = (DeuteriumDataValue) dataValue;
                value.rQ = nextGaussian(originalDataMap.get("rQ").get(i), value.rQError);
                value.rAP = nextGaussian(originalDataMap.get("rAP").get(i), value.rAPError);
            } else {
                throw new AssertionError("Unexpected subclass of `RelaxDataValue` detected!");
            }
        }
        // Force J-values to be computed
        data.jValues = null;
        return data;
    }

    public MolDataValues getOriginalData() {
        for (int i = 0; i < getNFields(); i++) {
            RelaxDataValue dataValue = data.dataValues.get(i);
            dataValue.R1 = originalDataMap.get("R1").get(i);
            dataValue.R2 = originalDataMap.get("R2").get(i);

            if (dataType.equals(R1R2NOEDataValue.class)) {
                R1R2NOEDataValue value = (R1R2NOEDataValue) dataValue;
                value.NOE = originalDataMap.get("NOE").get(i);
            } else if (dataType.equals(DeuteriumDataValue.class)) {
                DeuteriumDataValue value = (DeuteriumDataValue) dataValue;
                value.rQ = originalDataMap.get("rQ").get(i);
                value.rAP = originalDataMap.get("rAP").get(i);
            } else {
                throw new AssertionError("Unexpected subclass of `RelaxDataValue` detected!");
            }
        }
        // Force J-values to be computed
        data.jValues = null;
        return data;
    }
}
