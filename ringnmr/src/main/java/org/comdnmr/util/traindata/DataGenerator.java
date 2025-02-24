//ringnmr/src/main/java/org/comdnmr/util/traindata/DataGenerator.java
//Simon Hulse
//simonhulse@protonmail.com
//Last Edited: Tue 05 Nov 2024 03:47:56 PM EST

/**
 *
 *  @author simonhulse
 */

package org.comdnmr.util.traindata;

import org.comdnmr.util.DataUtil;
import org.comdnmr.util.traindata.samplers.*;
import org.comdnmr.util.traindata.dependants.*;

import org.comdnmr.eqnfit.CESTEquation;
import org.comdnmr.eqnfit.CPMGEquation;
import org.comdnmr.eqnfit.EquationType;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.text.StringSubstitutor;

import java.io.File;
import java.io.IOException;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class DataGenerator {

    static private final ObjectMapper objectMapper = new ObjectMapper();
    private final JsonNode config;
    private final int nDatasets;
    private final int nNoiseDuplicates;
    private final Pattern R1Matcher = Pattern.compile("^R1.*");
    private final Pattern R2Matcher = Pattern.compile("^R2.*");

    private Iterator<String[]> iterator;

    private Integer parSize;
    private Integer depSize;
    private Integer xSize;
    private Integer varSize;
    private List<Double> xValuesInterpolation;

    public String currentExperimentName;
    public String currentEnumName;
    String dataDirectoryTemplate;
    String dataPathTemplate;

    public static List<Double> getInterpolatedProfile(
        List<Double> yValues,
        List<Double> xValues,
        List<Double> xValuesInterpolation
    ) {
        int nValues = yValues.size();
        double[] xValuesArray = ArrayUtils.toPrimitive(xValues.toArray(new Double[nValues]));
        double[] yValuesArray = ArrayUtils.toPrimitive(yValues.toArray(new Double[nValues]));
        SplineInterpolator splineInterpolator = new SplineInterpolator();
        PolynomialSplineFunction spline = splineInterpolator.interpolate(xValuesArray, yValuesArray);

        List<Double> yInterpolated = new ArrayList<>();
        for (Double xValueInterpolation : xValuesInterpolation) {
            yInterpolated.add(spline.value(xValueInterpolation));
        }

        return yInterpolated;
    }

    public DataGenerator(String configFile) throws IOException, ReflectiveOperationException {
        config = objectMapper.readTree(new File(configFile));
        iterator = makeDataTypeIterator(config.get("data_types").fields());

        nDatasets = config.get("n_datasets").asInt();
        nNoiseDuplicates = config.get("n_noise_duplicates").asInt();
        dataDirectoryTemplate = config.get("enum_dir").asText();

        Map<String, String> sub = new HashMap<>();
        sub.put("root_dir", config.get("root_dir").asText());
        dataDirectoryTemplate = StringSubstitutor.replace(dataDirectoryTemplate, sub, "{", "}");

        dataPathTemplate = config.get("ring_data_path").asText();
    }

    public String getDataDirectory() {
        Map<String, String> sub = new HashMap<>();
        sub.put("experiment", currentExperimentName);
        sub.put("enum", currentEnumName);
        return StringSubstitutor.replace(dataDirectoryTemplate, sub, "{", "}");
    }

    public String getDataPath() {
        Map<String, String> sub = new HashMap<>();
        sub.put("enum_dir", getDataDirectory());
        return StringSubstitutor.replace(dataPathTemplate, sub, "{", "}");
    }

    private Iterator<String[]> makeDataTypeIterator(Iterator<Map.Entry<String, JsonNode>> dataTypes) {
        List<String[]> dataList = new ArrayList<>();

        while (dataTypes.hasNext()) {
            Map.Entry<String, JsonNode> info = dataTypes.next();
            String experimentName = info.getKey();
            Iterator<Map.Entry<String, JsonNode>> enums = info.getValue().get("enums").fields();
            while (enums.hasNext()) {
                String enumName = enums.next().getKey();
                String[] dataInfo = new String[2];
                dataInfo[0] = experimentName;
                dataInfo[1] = enumName;
                dataList.add(dataInfo);
            }
        }

        return dataList.iterator();
    }

    public boolean hasAnotherDataset() {
        return iterator.hasNext();
    }

    public void nextDataset() {
        if (iterator.hasNext()) {
            String[] info = iterator.next();
            currentExperimentName = info[0];
            currentEnumName = info[1];
        } else {
            currentExperimentName = null;
            currentEnumName = null;
        }
    }

    public Map<String, Object> generateDataset() {
        if (currentExperimentName == null) {
            return null;
        }

        // >>> CONSTRUCT CLASSES TO GENERATE EXPERIMENT AND MODEL PARAMETERS >>>
        JsonNode experimentInfo = config.get("data_types").get(currentExperimentName);
        // Sampler for the x-values (nuCPMG for CPMG, transmitter offset for CEST)
        Sampler xValuesSampler = getXValuesSampler(experimentInfo);
        // Sampler for the variable (B0 for CPMG, B1 field for CEST)
        Sampler variableSampler = getVariableSampler(experimentInfo);
        // Constants are experiment parameters that are fixed for a given set
        // of profiles. For example, tau in CPMG.
        List<Sampler> constantSamplers = getConstantSamplers(experimentInfo);
        // Generators for experiment parameters that are needed, but are
        // directly related to the varaible.
        //
        // I included this specifically because `CPMGMQ` requires both the 1H
        // field and the heteronucleus field, and the heteronucleus field is of
        // course related to the 1H field by a scaling factor based on the
        // relevant gyromagnetic ratios.
        List<DependantGenerator> dependantGenerators = getDependantGenerators(experimentInfo);
        // Samplers for the model parameters, for example R2, kEx, pA, deltaCPPM in CPMGSLOW
        List<Sampler> parameterSamplers = getParameterSamplers(experimentInfo);
        // x-values used to interpolate the data
        xValuesInterpolation = getXValuesInterpolation(experimentInfo);
        // <<<

        setVarSize(variableSampler);
        setParSize(parameterSamplers);
        setDepSize(dependantGenerators);
        setXSize(constantSamplers);

        // >>> CONSTRUCT LISTS TO STORE DATASETS >>>
        DataList<List<Double>> xValueList = new DataList<>(xValuesSampler.name);

        DataList<List<List<Double>>> profileList = new DataList<>("profile");
        DataList<List<List<Double>>> profileNoisyList = new DataList<>("profile_noisy");
        DataList<List<List<Double>>> profileInterpolatedList = new DataList<>("profile_interpolation");

        List<DataList<Double>> parameterLists = initializeListsFromSamplers(parameterSamplers);
        DataList<List<Double>> variableList = new DataList<>(variableSampler.name);
        List<DataList<List<Double>>> dependantLists = initializeListsFromGenerators(dependantGenerators);
        List<DataList<Double>> constantLists = initializeListsFromSamplers(constantSamplers);
        // <<<

        EquationType cls = fetchEquationClass();

        int[][] map = makeMaps(parameterSamplers);

        for (int n = 0; n < this.nDatasets; n++) {
            List<Double> xValues = xValuesSampler.sample();
            List<Double> variables = variableSampler.sample();
            List<Double> constants = sampleConstants(constantSamplers);
            List<List<Double>> dependants = generateDependants(variables, dependantGenerators);

            double[] par = makePar(parameterSamplers, variables);
            double[][][] x = makeX(xValues, variables, constants, dependants);
            List<List<Double>> profiles = makeProfiles(cls, par, map, x);

            for (int p = 0; p < nNoiseDuplicates; p++) {
                List<List<Double>> noise = makeNoise(experimentInfo, profiles);
                List<List<Double>> profilesNoisy = sumTwoDLists(profiles, noise);
                List<List<Double>> profilesInterpolated = makeProfilesInterpolated(profilesNoisy, xValues);

                profileList.add(profiles);
                profileNoisyList.add(profilesNoisy);
                profileInterpolatedList.add(profilesInterpolated);

                xValueList.add(xValues);
                variableList.add(variables);

                appendToListOfDataLists(dependantLists, dependants);
                appendToListOfDataLists(parameterLists, par);
                appendToListOfDataLists(constantLists, constants);
                printMilestoneMessage(n, p);
            }
        }

        printCompleteMessage();

        return bundleData(
            parameterLists,
            constantLists,
            variableList,
            dependantLists,
            xValueList,
            profileList,
            profileNoisyList,
            profileInterpolatedList);
    }

    void printMilestoneMessage(int n, int p) {
        int datasetNumber = n * nNoiseDuplicates + p;
        if ((datasetNumber + 1) % 10000 == 0 && (datasetNumber + 1) != nDatasets * nNoiseDuplicates) {
            System.out.println(
                String.format(
                    "[%s - %s]: %d/%d datasets created",
                    currentExperimentName,
                    currentEnumName,
                    datasetNumber + 1,
                    nDatasets * nNoiseDuplicates
                )
            );
        }
    }

    void printCompleteMessage() {
        System.out.println(
            String.format(
                "[%s - %s]: Finished making %d datasets",
                currentExperimentName,
                currentEnumName,
                nDatasets * nNoiseDuplicates
            )
        );
    }

    private Sampler getXValuesSampler(JsonNode experimentInfo) {
        return fetchSampler(experimentInfo.get("x_values"));
    }

    private Sampler getVariableSampler(JsonNode experimentInfo) {
        return fetchSampler(experimentInfo.get("variable"));
    }

    private List<DependantGenerator> getDependantGenerators(JsonNode experimentInfo) {
        JsonNode dependantInfo = experimentInfo.get("variable").get("dependants");
        if (dependantInfo == null) {
            return null;
        }

        Iterator<JsonNode> dependantIterator = dependantInfo.elements();
        List<DependantGenerator> dependantGenerators = new ArrayList<DependantGenerator>();

        while (dependantIterator.hasNext()) {
            dependantGenerators.add(fetchDependantGenerator(dependantIterator.next()));
        }

        return dependantGenerators;
    }

    private DependantGenerator fetchDependantGenerator(JsonNode dependantInfo) {
        DependantGenerator generator = null;

        String name = dependantInfo.get("name").asText();
        String generatorType = dependantInfo.get("relation").asText();

        switch (generatorType) {
            case "Multiply":
                Sampler sampler = fetchSampler(dependantInfo);
                generator = new Multiplier(name, sampler);
                break;

            default: throw new RuntimeException("Unknown dependant specification");
        }

        return generator;
    }

    private List<Sampler> getConstantSamplers(JsonNode experimentInfo) {
        List<JsonNode> constantInfoList = objectMapper.convertValue(
            experimentInfo.get("constants"),
            new TypeReference<ArrayList<JsonNode>>() { }
        );

        List<Sampler> constantSamplers = new ArrayList<>();
        for (JsonNode constantInfo : constantInfoList) {
            constantSamplers.add(fetchSampler(constantInfo));
        }
        return constantSamplers;
    }

    private List<Sampler> getParameterSamplers(JsonNode experimentInfo) {
        Iterator<JsonNode> parameterIterator = experimentInfo.get("enums").get(currentEnumName).elements();
        List<Sampler> parameterSamplers = new ArrayList<>();

        while (parameterIterator.hasNext()) {
            JsonNode parameterInfo = parameterIterator.next();
            String parameterName = parameterInfo.get(0).asText();
            double parameterMin = parameterInfo.get(1).asDouble();
            double parameterMax = parameterInfo.get(2).asDouble();
            parameterSamplers.add(new RandomDoubleSampler(parameterName, parameterMin, parameterMax, 1, 1));
        }

        return parameterSamplers;
    }

    private void setVarSize(Sampler variableSampler) {
        varSize = variableSampler.sample().size();
    }

    private int isRelaxationRateSampler(Sampler parameterSampler) {
        if (R1Matcher.matcher(parameterSampler.name).matches()) {
            return 1;
        } else if (R2Matcher.matcher(parameterSampler.name).matches()) {
            return 2;
        } else {
            return 0;
        }
    }

    private void setParSize(List<Sampler> parameterSamplers) {
        parSize = 0;
        for (Sampler parameterSampler : parameterSamplers) {
            // If the sampler is for a relaxation rate we need one parameter per profile.
            // Otherwise, the parameter is global across profiles.
            parSize += (isRelaxationRateSampler(parameterSampler) > 0) ? varSize : 1;
        }
    }

    private void setDepSize(List<DependantGenerator> dependantGenerators) {
        depSize = dependantGenerators != null ? dependantGenerators.size() : 0;
    }

    private void setXSize(List<Sampler> constantSamplers) {
        // The x-values provided per sample are:
        // 0 --> x-value
        // 1, ... --> dependencies
        // 1 + n(dependencies) --> variable
        // 2 + n(dependencies), ... --> constants
        xSize = 2 + depSize + constantSamplers.size();
    }

    private List<Double> getXValuesInterpolation(JsonNode experimentInfo) {
        return fetchSampler(experimentInfo.get("x_values_interpolation")).sample();
    }

    private List<DataList<Double>> initializeListsFromSamplers(List<Sampler> samplers) {
       List<DataList<Double>> lists = new ArrayList<>();
        for (Sampler sampler : samplers) {
            String name = sampler.getName();
            if (isRelaxationRateSampler(sampler) > 0) {
                String nameTemplate = "%s<%d>";
                for (int i = 1; i <= varSize; i++) {
                    lists.add(new DataList<Double>(String.format(nameTemplate, name, i)));
                }
            } else {
                lists.add(new DataList<Double>(name));
            }
        }
        return lists;
    }

    private List<DataList<List<Double>>> initializeListsFromGenerators(List<DependantGenerator> generators) {
       List<DataList<List<Double>>> lists = new ArrayList<>();
        for (DependantGenerator generator : generators) {
            String name = generator.getName();
            lists.add(new DataList<List<Double>>(name));
        }
        return lists;
    }

    private List<Double> sampleConstants(List<Sampler> constantSamplers) {
        List<Double> constants = new ArrayList<>();
        for (Sampler constantSampler : constantSamplers) {
            constants.add(constantSampler.sample().get(0));
        }
        return constants;
    }

    private List<List<Double>> generateDependants(List<Double> variables, List<DependantGenerator> dependantGenerators) {
        if (dependantGenerators == null) {
            return null;
        }
        List<List<Double>> dependants = new ArrayList<>();
        for (DependantGenerator dependantGenerator : dependantGenerators) {
            dependantGenerator.updateState();
            List<Double> dependantList = new ArrayList<>();
            for (double variable : variables) {
                dependantList.add(dependantGenerator.fetch(variable));
            }
            dependants.add(dependantList);
        }

        return dependants;
    }

    private int[][] makeMaps(List<Sampler> parameterSamplers) {
        int[][] maps = new int[varSize][parameterSamplers.size()];
        List<int[]> indices = new ArrayList<>();
        int idx = 0;
        for (Sampler parameterSampler : parameterSamplers) {
            int[] array;
            if (isRelaxationRateSampler(parameterSampler) > 0) {
                array = IntStream.range(idx, idx + varSize).toArray();
                idx += varSize;
            } else {
                array = new int[]{idx++};
            }
            indices.add(array);
        }
        for (int i = 0; i < varSize; i++) {
            for (int j = 0; j < indices.size(); j++) {
                int[] idxs = indices.get(j);
                maps[i][j] = idxs[i % idxs.length];
            }
        }

        return maps;
    }

    private double[] makePar(List<Sampler> parameterSamplers, List<Double> variables) {
        double[] par = new double[parSize];
        int idx = 0;
        double kEx = 0.0;
        double pA = 0.0;
        double thold = 0.0;
        double maxVariable = Collections.max(variables);
        for (Sampler parameterSampler : parameterSamplers) {
            double parameter = parameterSampler.sample().get(0);

            // TODO: more elegant solution to this?
            if (idx == 0) {
                kEx = parameter;
            }
            else if (idx == 1) {
                pA = parameter;
            }
            else {
                switch (currentEnumName) {
                    case "CPMGSLOW":
                        if (idx == 2 + variables.size()) {
                            // Ensure slow regime is adhered to
                            // Only accept deltaXPPM when it is greater than:
                            // 5 kEx / 2 pi B0max
                            thold = (5.0 * kEx) / (2.0 * Math.PI * maxVariable);
                            while (parameter < thold) {
                                parameter = parameterSampler.sample().get(0);
                            }
                        }
                        break;

                    case "CPMGFAST":
                        if (idx == 1 + variables.size()) {
                            // Ensure fast regime is adhered to
                            // Only accept deltaPPMmin when it is less than:
                            // 5 kEx / 2 pi B0max
                            // N.B. Factor of 1/2 is for max(sqrt(pa pb))
                            thold = (5.0 * kEx) / (2.0 * Math.PI * maxVariable);
                            while (parameter > thold) {
                                parameter = parameterSampler.sample().get(0);
                            }
                        }
                        break;
                }
            }

            int parameterSamplerType = isRelaxationRateSampler(parameterSampler);
            if (parameterSamplerType == 0) {
                par[idx++] = parameter;
            }
            else {
                double rateScale = (parameterSamplerType == 1)
                                 ? (config.get("R1_scale_factor").asDouble())
                                 : (config.get("R2_scale_factor").asDouble());
                List<Double> rates = RelaxationRateGenerator.generate(parameter, variables, rateScale);
                for (double rate : rates) {
                    par[idx++] = rate;
                }
            }
        }

        return par;
    }

    private double[][][] makeX(
        List<Double> xValues,
        List<Double> variables,
        List<Double> constants,
        List<List<Double>> dependants
    ) {
        int xValSize = xValues.size();
        double[][][] x = new double[varSize][xValSize][xSize];

        for (int varIdx = 0; varIdx < varSize; varIdx++) {
            double variable = variables.get(varIdx);
            for (int xValIdx = 0; xValIdx < xValSize; xValIdx++) {
                int idx = 0;
                x[varIdx][xValIdx][idx++] = xValues.get(xValIdx);
                if (dependants != null) {
                    for (List<Double> dependantList : dependants) {
                        double dep = dependantList.get(varIdx);
                        x[varIdx][xValIdx][idx++] = dep;
                    }
                }
                x[varIdx][xValIdx][idx++] = variable;
                for (double constant : constants) {
                    x[varIdx][xValIdx][idx++] = constant;
                }
            }
        }

        return x;
    }

    private void appendToListOfDataLists(List<DataList<Double>> lists, double[] data) {
        for (int i = 0; i < data.length; i++) {
            lists.get(i).add(data[i]);
        }
    }

    private <T> void appendToListOfDataLists(List<DataList<T>> lists, List<T> data) {
        for (int i = 0; i < lists.size(); i++) {
            lists.get(i).add(data.get(i));
        }
    }

    private List<List<Double>> makeProfiles(EquationType cls, double[] par, int[][] map, double[][][] x) {
        List<List<Double>> profiles = new ArrayList<>();
        for (int varIdx = 0; varIdx < varSize; varIdx++) {
            List<Double> profile = new ArrayList<>();
            for (int xValueIdx = 0; xValueIdx < x[0].length; xValueIdx++) {
                profile.add(cls.calculate(par, map[varIdx], x[varIdx][xValueIdx], 0));
            }
            profiles.add(profile);
        }
        return profiles;
    }

    private List<List<Double>> makeProfilesInterpolated(List<List<Double>> profiles, List<Double> xValues) {
        if (xValues.equals(xValuesInterpolation)) {
            return profiles;
        }
        List<List<Double>> profilesInterpolated = new ArrayList<>();
        for (List<Double> profile : profiles) {
            profilesInterpolated.add(DataGenerator.getInterpolatedProfile(profile, xValues, xValuesInterpolation));
        }
        return profilesInterpolated;
    }

    private Map<String, Object> bundleData(
        List<DataList<Double>> parameterLists,
        List<DataList<Double>> constantLists,
        DataList<List<Double>> variableList,
        List<DataList<List<Double>>> dependantLists,
        DataList<List<Double>> xValueList,
        DataList<List<List<Double>>> profileList,
        DataList<List<List<Double>>> profileNoisyList,
        DataList<List<List<Double>>> profileInterpolationList
    ) {
        Map<String, Object> data = new HashMap<>();

        List<Map<String, Object>> parameterMapList = new ArrayList<>();
        for (DataList<Double> parameterList : parameterLists) {
            parameterMapList.add(parameterList.asHashMap());
        }

        List<Map<String, Object>> constantMapList = new ArrayList<>();
        for (DataList<Double> constantList : constantLists) {
            constantMapList.add(constantList.asHashMap());
        }

        List<Map<String, Object>> dependantMapList = new ArrayList<>();
        for (DataList<List<Double>> dependantList : dependantLists) {
            dependantMapList.add(dependantList.asHashMap());
        }

        data.put("parameters", parameterMapList);
        data.put("constants", constantMapList);
        data.put("dependants", dependantMapList);
        data.put("variable", variableList.asHashMap());
        data.put("x_values", xValueList.asHashMap());
        data.put("profile", profileList.asHashMap());
        data.put("profile_noisy", profileNoisyList.asHashMap());
        data.put("profile_interpolation", profileInterpolationList.asHashMap());

        return data;
    }


    private Sampler fetchSampler(JsonNode info) {
        Sampler sampler = null;

        String name = info.get("name").asText();
        JsonNode samplerInfo = info.get("sampling");

        switch (samplerInfo.get("type").asText()) {
            case "RandomDouble":
                sampler = new RandomDoubleSampler(
                    name,
                    samplerInfo.get("minValue").asDouble(),
                    samplerInfo.get("maxValue").asDouble(),
                    samplerInfo.get("minSamples").asInt(),
                    samplerInfo.get("maxSamples").asInt()
                );
                break;

            case "Choice":
                sampler = new RandomChoiceSampler(
                    name,
                    objectMapper.convertValue(samplerInfo.get("options"), ArrayList.class),
                    samplerInfo.get("minSamples").asInt(),
                    samplerInfo.get("maxSamples").asInt()
                );
                break;

            case "RandomLinspace":
                sampler = new RandomLinspaceSampler(
                    name,
                    samplerInfo.get("minFirst").asDouble(),
                    samplerInfo.get("maxFirst").asDouble(),
                    samplerInfo.get("minLast").asDouble(),
                    samplerInfo.get("maxLast").asDouble(),
                    samplerInfo.get("minSamples").asInt(),
                    samplerInfo.get("maxSamples").asInt()
                );
                break;


            case "RandomLogspace":
                sampler = new RandomLogspaceSampler(
                    name,
                    samplerInfo.get("minFirst").asDouble(),
                    samplerInfo.get("maxFirst").asDouble(),
                    samplerInfo.get("scale").asDouble(),
                    samplerInfo.get("minSamples").asInt(),
                    samplerInfo.get("maxSamples").asInt()
                );
                break;

            case "RandomSegmentedLinspace":
                sampler = new RandomSegmentedLinspaceSampler(
                    name,
                    objectMapper.convertValue(samplerInfo.get("minFirsts"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("maxFirsts"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("minLasts"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("maxLasts"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("minSampless"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("maxSampless"), ArrayList.class)
                );
                break;

            case "Explicit":
                sampler = new ExplicitSampler(
                    name,
                    objectMapper.convertValue(samplerInfo.get("value"), ArrayList.class)
                );
                break;

            default: throw new RuntimeException("Unrecognised sampler.");
        }

        return sampler;
    }

    private EquationType fetchEquationClass() {
        EquationType cls = null;

        switch (currentExperimentName) {
            case "CPMGEquation":
                switch (currentEnumName) {
                    case "CPMGMQ": cls = CPMGEquation.CPMGMQ; break;
                    case "CPMGFAST": cls = CPMGEquation.CPMGFAST; break;
                    case "CPMGSLOW": cls = CPMGEquation.CPMGSLOW; break;
                }
                break;

            case "CESTEquation":
                switch (currentEnumName) {
                    case "EXACT0": cls = CESTEquation.EXACT0; break;
                }
                break;
        }

        if (cls == null) {
            throw new RuntimeException("Unrecognised model.");
        }

        return cls;
    }

    private List<List<Double>> makeNoise(JsonNode experimentInfo, List<List<Double>> profiles) {
        double maxVariance = experimentInfo.get("variance_max").asDouble();
        int nValues = profiles.get(0).size();
        RandomDoubleSampler varianceSampler = new RandomDoubleSampler(0.0, maxVariance, 1, 1);
        List<List<Double>> noises = new ArrayList<>();
        for (int i = 0; i < varSize; i++) {
            noises.add(sampleAWGN(varianceSampler.sampleUnwrapped(), nValues));
        }
        return noises;
    }

    private List<Double> sampleAWGN(double variance, int nValues) {
        NormalDistribution sampler = new NormalDistribution(0.0, Math.sqrt(variance));
        List<Double> noise = new ArrayList<>();
        for (int i = 0; i < nValues; i++) {
            noise.add(sampler.sample());
        }
        return noise;
    }

    private List<List<Double>> sumTwoDLists(List<List<Double>> l1, List<List<Double>> l2) {
        List<List<Double>> sum = new ArrayList<>();
        for (Map.Entry<List<Double>, List<Double>> vecPair : zipLists(l1, l2)) {
            List<Double> vec1 = vecPair.getKey();
            List<Double> vec2 = vecPair.getValue();
            List<Double> vecSum = new ArrayList<>();
            for (Map.Entry<Double, Double> valPair : zipLists(vec1, vec2)) {
                double val1 = valPair.getKey();
                double val2 = valPair.getValue();
                vecSum.add(val1 + val2);
            }
            sum.add(vecSum);
        }
        return sum;
    }

    private static <A, B> List<Map.Entry<A, B>> zipLists(List<A> l1, List<B> l2) {
        return IntStream.range(0, Math.min(l1.size(), l2.size()))
            .mapToObj(i -> Map.entry(l1.get(i), l2.get(i)))
            .collect(Collectors.toList());
    }

    private Integer interpolationSize() {
        if (currentExperimentName == null) {
            return null;
        }
        return xValuesInterpolation.size();
    }

}

class RelaxationRateGenerator {

    static List<Double> generate(double baseRate, List<Double> variable, double maxMultiplier) {
        double minVariable = variable.get(0);
        double currMax = baseRate;
        RandomDoubleSampler sampler = new RandomDoubleSampler(0.0, maxMultiplier, 1, 1);
        List<Double> rates = new ArrayList<>();
        rates.add(baseRate);
        int i = 1;
        while (i < variable.size()) {
            double scale = (variable.get(i) - minVariable) / minVariable;
            double randomMultiplier = sampler.sampleUnwrapped();

            // Ensure that rate for higher fields is always larger.
            // This is at least the case for R2... may have to rethink this for R1.
            double candidate = baseRate + randomMultiplier * scale;
            if (candidate > currMax) {
                currMax = candidate;
                rates.add(candidate);
                i++;
            }
        }
        return rates;
    }
}

