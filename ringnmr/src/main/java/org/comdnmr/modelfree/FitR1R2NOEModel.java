/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;
import org.nmrfx.chemistry.*;
import org.nmrfx.chemistry.relax.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Class for fitting the model-free formalism to spectral density data
 * generated using reduced spectral density mapping.
 *
 * @author brucejohnson
 * @author simonhulse
 */
public class FitR1R2NOEModel extends FitModel {
    Map<String, MolDataValues> molData = null;

    private static final String[] PARAMETER_NAMES = {"Tau_e", "Sf2", "Tau_f", "Ss2", "Tau_s"};

    public int getNPars() { return PARAMETER_NAMES.length; }

    public void setData(Map<String, MolDataValues> molData) {
        this.molData = molData;
    }

    Map<String, MolDataValues> getData(boolean requireCoords) {
        Map<String, MolDataValues> molDataValues = new HashMap<>();
        MoleculeBase mol = MoleculeFactory.getActive();
        if (mol != null) {
            for (Polymer polymer : mol.getPolymers()) {
                for (Residue residue : polymer.getResidues()) {
                    var hAtom = residue.getAtom("H");
                    var nAtom = residue.getAtom("N");
                    if (nAtom == null) {
                        continue;
                    }
                    Vector3D vec = null;
                    if (requireCoords) {
                        if (hAtom != null) {
                            var hPoint = hAtom.getPoint();
                            var nPoint = nAtom.getPoint();
                            if ((hPoint == null) || (nPoint == null)) {
                                continue;
                            } else {
                                vec = hPoint.subtract(nPoint);
                            }
                        } else {
                            continue;
                        }
                    }
                    var relaxData = nAtom.getRelaxationData();

                    if (relaxData != null) {
                        MolDataValues molData = vec == null
                                ? new MolDataValues(nAtom)
                                : new MolDataValues(nAtom, vec.toArray());
                        Set<Long> fields = new HashSet<>();
                        for (var entry : relaxData.entrySet()) {
                            RelaxationData data = entry.getValue();
                            RelaxationSet relaxationSet = entry.getKey();
                            if (relaxationSet.active()) {
                                fields.add(Math.round(relaxationSet.field()));
                            }
                        }
                        for (var field : fields) {
                            Double r1 = null;
                            Double r1Error = null;

                            Double r2 = null;
                            Double r2Error = null;

                            Double noe = null;
                            Double noeError = null;

                            for (var entry : relaxData.entrySet()) {
                                RelaxationData data = entry.getValue();
                                if (data.getResonanceSource().deleted()) {
                                    break;
                                }
                                RelaxationSet relaxationSet = entry.getKey();
                                if (!relaxationSet.active()) {
                                    continue;
                                }
                                if (Math.round(entry.getKey().field()) == field) {
                                    switch (entry.getKey().relaxType()) {
                                        case R1 -> {
                                            r1 = data.getValue();
                                            r1Error = data.getError();
                                        }
                                        case R2 -> {
                                            r2 = data.getValue();
                                            r2Error = data.getError();
                                        }
                                        case NOE -> {
                                            noe = data.getValue();
                                            noeError = data.getError();
                                        }
                                    }
                                }
                            }
                            if ((r1 != null) && (r2 != null) && (noe != null)) {
                                RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, "H", "N");

                                R1R2NOEDataValue dValue = new R1R2NOEDataValue(molData, r1, r1Error, r2, r2Error, noe, noeError, relaxObj);
                                molData.addData(dValue);
                            }
                        }
                        if (!molData.dataValues.isEmpty()) {
                            molDataValues.put(molData.specifier, molData);
                        }
                    }
                }
            }
        }
        return molDataValues;
    }

    public  Map<String, ModelFitResult> testIsoModel() {
        if ((molData == null) || (molData.isEmpty())) {
            molData = getData(false);
        }
        if (searchKey != null) {
            if (molData.containsKey(searchKey)) {
                var keepVal = molData.get(searchKey);
                molData.clear();
                molData.put(searchKey, keepVal);
            }
        }

        if (!molData.isEmpty()) {
            if (fitSpec.tauMNeedsComputing) fitSpec.setTauM(estimateTau(molData).get("tau"));

            Map<String, ModelFitResult> results = new HashMap<>();
            molData
                .entrySet()
                .stream()
                .sorted(Comparator.comparing(Map.Entry::getKey))
                // .parallel()
                .forEach(
                    residue ->
                        results.put(
                            residue.getKey(),
                            fitSpec.fit(residue.getKey(), residue.getValue())
                        )
                );
            return results;
        } else {
            throw new IllegalStateException("No relaxation data to analyze. Need T1, T2 and NOE datasets");
        }

    }

    public Map<String, Double> estimateTau() {
        Map<String, MolDataValues> molData = getData(false);
        return estimateTau(molData);
    }

    Map<String, Double> estimateTau(Map<String, MolDataValues> molData) {
        Map<Long, List<RelaxDataValue>> map = new HashMap<>();
        for (var entry : molData.entrySet()) {
            MolDataValues value = entry.getValue();
            for (RelaxDataValue rlxValue : value.getData()) {
                long sfMHz = Math.round(rlxValue.relaxObj.getSF() / 1.0e6);
                List<RelaxDataValue> values = map.computeIfAbsent(sfMHz, k -> new ArrayList<>());
                values.add(rlxValue);
            }
        }
        int max = 0;
        long maxMHz = 0;
        for (long sfMHz : map.keySet()) {
            int size = map.get(sfMHz).size();
            if (size > max) {
                maxMHz = sfMHz;
                max = size;
            }
        }
        if (max > 0) {
            return CorrelationTime.estimateTau(maxMHz, "N", map.get(maxMHz));
        } else {
            return Collections.emptyMap();
        }
    }


    public  Map<String, ModelFitResult> testModels(Map<String, MolDataValues> molData, List<String> modelNames) {
        Random random = new Random();
        if (tau == null) tau = estimateTau(molData).get("tau");
        AtomicInteger counts = new AtomicInteger();
        int n = molData.entrySet().size();
        Map<String, ModelFitResult> results = new ConcurrentHashMap<>();
        MoleculeBase moleculeBase = MoleculeFactory.getActive();
        Map<String, OrderParSet> orderParSetMap = moleculeBase.orderParSetMap();

        for (var modelName : modelNames) {
            OrderParSet orderParSet = orderParSetMap.computeIfAbsent("order_parameter_list_" + modelName, k -> new OrderParSet(k));
        }
        OrderParSet orderParSet = orderParSetMap.computeIfAbsent("order_parameter_list_best", k -> new OrderParSet(k));

        molData.entrySet().stream().sorted(Comparator.comparing(Map.Entry::getKey)).parallel().forEach(e -> {

            updateProgress((double) counts.get() / n);
            if (cancelled.get()) {
                return;
            }
            MolDataValues resData = e.getValue();
            String key = e.getKey();
            if (!resData.getData().isEmpty()) {
                Optional<ModelFitResult> result;
                switch (bootstrapMode) {
                    case PARAMETRIC:
                        result = testModels(
                            orderParSetMap, resData, key, modelNames, random);
                        break;
                    case AGGREGATE:
                    case BAYESIAN:
                        result = testModelsWithBootstrapAggregation(
                            orderParSetMap, resData, key, modelNames, random);
                        break;
                    default: throw new AssertionError("Unreachable");
                }
                result.ifPresent(o -> results.put(key, o));
            }
            int iCount = counts.incrementAndGet();
        });
        return results;
    }

    private OrderPar makeOrderPar(
        OrderParSet orderParSet,
        MolDataValues resData,
        String key,
        Score bestScore,
        MFModelIso bestModel,
        double[][] repData
    ) {
        ResonanceSource resSource = new ResonanceSource(resData.atom);
        Atom atom = resData.atom;
        var parNames = bestModel.getParNames();
        double[] pars = bestScore.getPars();

        OrderPar orderPar = new OrderPar(orderParSet, resSource, bestScore.rss, bestScore.nValues, bestScore.nPars, bestModel.getName());
        for (int iPar = 0; iPar < parNames.size(); iPar++) {
            String parName = parNames.get(iPar);
            double parValue = pars[iPar];
            Double parError = null;
            if (repData != null) {
                DescriptiveStatistics sumStat = new DescriptiveStatistics(repData[iPar]);
                parError = sumStat.getStandardDeviation();
            }
            orderPar = orderPar.set(parName, parValue, parError);
        }
        if (!bestModel.fitTau()) {
            orderPar = orderPar.set("Tau_e", bestModel.getTau(), 0.0);
        }

        orderPar = orderPar.set("model", (double) bestModel.getNumber(), null);
        if ((bestModel instanceof MFModelIso2sf) && useLambda()) {
            orderPar = orderPar.setModel();
        }
        atom.addOrderPar(orderParSet, orderPar);
        double[][] jValues = resData.getJValues();
        SpectralDensity spectralDensity = new SpectralDensity(key, jValues);
        atom.addSpectralDensity(key, spectralDensity);
        return orderPar;

    }

    public Optional<ModelFitResult> testModels(
        Map<String, OrderParSet> orderParSetMap,
        MolDataValues resData,
        String key,
        List<String> modelNames,
        Random random
    ) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        molDataRes.put(key, resData);

        boolean localFitTau = fitTau(molDataRes);
        double localTauFraction = (localFitTau) ? tauFraction : 0.0;

        MFModelIso bestModel = null;
        Score bestScore = null;
        double lowestAIC = Double.MAX_VALUE;

        for (var modelName : modelNames) {
            MFModelIso model = MFModelIso.buildModel(
                modelName, localFitTau, tau, localTauFraction, fitExchange);
            resData.setTestModel(model);
            Score score = tryModel(molDataRes, model, localTauFraction, localFitTau, random);
            OrderParSet orderParSet = orderParSetMap.get("order_parameter_list_"+ modelName);
            double[][] repData = null;
            if (nReplicates > 2) {
                double[] pars = score.getPars();
                repData = replicates(molDataRes, model, localTauFraction, localFitTau, pars, random);
                OrderPar orderPar = makeOrderPar(orderParSet, resData, key, score, model, repData);
            }
            if (useLambda || (modelNames.size() == 1) || (score.aicc().isPresent() && ((bestScore == null) || (score.aicc().get() < lowestAIC))))  {
                lowestAIC = score.aicc().isPresent() ? score.aicc().get() : null;
                bestModel = model;
                bestScore = score;
            }

        }
        Optional<ModelFitResult> result = Optional.empty();
        double[][] replicateData;
        if (bestScore != null) {
            double[] pars = bestScore.getPars();
            var parNames = bestModel.getParNames();
            replicateData = null;
            if (nReplicates > 2) {
                replicateData = replicates(molDataRes, bestModel, localTauFraction, localFitTau, pars, random);
            }
            OrderParSet orderParSet = orderParSetMap.get("order_parameter_list_best");
            OrderPar orderPar = makeOrderPar(orderParSet, resData, key, bestScore, bestModel, replicateData);
            ModelFitResult modelFitResult = new ModelFitResult(
                orderPar,
                replicateData,
                null,
                new Score[]{bestScore}
            );
            result = Optional.of(modelFitResult);
        }
        return result;
    }

    /**
     * Perform model-free fitting using <a href="https://mr.copernicus.org/articles/2/251/2021/">bootstrap aggregation</a>.
     *
     * @param orderParSetMap a map containing order parameter sets indexed by string keys.
     * @param resData instance of {@link MolDataValues} containing experimental data for fitting.
     * @param key a string key representing the identifier for spectral density data.
     * @param modelNames a list of model names to be tested during the fitting. Allowed elements are `1`, `1f`, `1s`, `2s`, `2sf`.
     * @param random an instance of {@link Random} for generating random numbers during sampling.
     * @return an {@link Optional<ModelFitResult>} containing the result of the model fitting, or empty if cancelled or no valid result is found.
     */
    public Optional<ModelFitResult> testModelsWithBootstrapAggregation(
        Map<String, OrderParSet> orderParSetMap,
        MolDataValues resData,
        String key,
        List<String> modelNames,
        Random random
    ) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        molDataRes.put(key, resData);

        boolean localFitTau = fitTau(molDataRes);
        double localTauFraction = (localFitTau) ? tauFraction : 0.0;

        BootstrapSampler sampler;
        switch (bootstrapMode) {
            case AGGREGATE:
                sampler = new NonparametricSampler(resData);
                break;
            case BAYESIAN:
                sampler = new BayesianSampler(resData);
                break;
            default:
                throw new AssertionError("Unreachable");
        }

        Score[] scores = new Score[nReplicates];
        MFModelIso[] models = new MFModelIso[nReplicates];
        // Using the same index variables as Crawley and Palmer
        // i: boostrap repliacte
        // j: spectral density datapoint
        // k: parameter index
        for (int i = 0; i < nReplicates; i++) {
            if (cancelled.get()) {
                return Optional.empty();
            }

            MolDataValues replicateData = sampler.sample();
            molData = new HashMap<>() {{ put(key, replicateData); }};

            // Fit each model and retain the one which has the lowest AICc
            Optional<Pair<MFModelIso, Score>> modelScore = Optional.empty();
            for (var modelName : modelNames) {
                MFModelIso testModel = MFModelIso.buildModel(
                        modelName, localFitTau, tau, localTauFraction, fitExchange);
                resData.setTestModel(testModel);
                Score testScore = tryModel(molDataRes, testModel, localTauFraction, localFitTau, random);
                if (modelScore.isEmpty() || testScore.aicc().get() < modelScore.get().getValue().aicc().get()) {
                    modelScore = Optional.of(Pair.of(testModel, testScore));
                }
            }

            if (modelScore.isEmpty()) return Optional.empty();
            models[i] = modelScore.get().getKey();
            scores[i] = modelScore.get().getValue();
        }

        MFModelIso model2sf = MFModelIso.buildModel("2sf", localFitTau, tau, localTauFraction, fitExchange);

        ResonanceSource resSource = new ResonanceSource(resData.atom);
        double rssSum = 0.0;
        for (int i = 0; i < nReplicates; i++) {
            rssSum += scores[i].rss;
        }
        double rss = rssSum /= nReplicates;
        OrderParSet orderParSet = orderParSetMap.get("order_parameter_list_best");
        OrderPar orderPar = new OrderPar(
            orderParSet,
            resSource,
            rss,
            scores[0].nValues,
            getNPars(),
            model2sf.getName());
        orderPar = orderPar.set("model", (double) model2sf.getNumber(), null);

        computeBaggingStats(orderPar, models, scores);

        Atom atom = resData.atom;
        atom.addOrderPar(orderParSet, orderPar);
        double[][] jValues = resData.getJValues();
        SpectralDensity spectralDensity = new SpectralDensity(key, jValues);
        atom.addSpectralDensity(key, spectralDensity);
        resData.setTestModel(model2sf);

        ModelFitResult result = new ModelFitResult(orderPar, makeBaggingParameters(models, scores), null, scores);
        return Optional.of(result);
    }

    private double[][] makeBaggingParameters(MFModelIso[] models, Score[] scores) {
        int nReplicates = scores.length;
        double[][] parameters = new double[getNPars()][nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            double[] replicateParameters = models[i].getStandardPars(scores[i].getPars());
            for (int k = 0; k < getNPars(); k++) {
                parameters[k][i] = replicateParameters[k];
            }
        }
        return parameters;
    }

    private double[] makeBaggingWeightMeans(Score[] scores) {
        int nReplicates = scores.length;
        int nSpecdens = scores[0].getWeights().length;
        double[] weightMeans = new double[nSpecdens];
        for (int i = 0; i < nReplicates; i++) {
            double[] replicateWeights = scores[i].getWeights();
            for (int j = 0; j < nSpecdens; j++) {
                weightMeans[j] += replicateWeights[j];
            }
        }
        for (int j = 0; j < nSpecdens; j++) {
            weightMeans[j] /= nReplicates;
        }
        return weightMeans;
    }

    private void computeBaggingStats(OrderPar orderPar, MFModelIso[] models, Score[] scores) {
        double[] weightMeans = makeBaggingWeightMeans(scores);
        double[][] parameters = makeBaggingParameters(models, scores);

        int nSpecdens = scores[0].getWeights().length;
        int nReplicates = scores.length;
        for (int k = 0; k < getNPars(); k++) {
            String name = PARAMETER_NAMES[k];

            // Get parameter means
            double mean = 0.0;
            for (int i = 0; i < nReplicates; i++) {
                mean += parameters[k][i];
            }
            mean /= nReplicates;

            // Get smoothed parameter errors
            // Eqs 17-19 in Crawley and Palmer's paper
            // N.B. there is an erroneous 1/N in eq 19 which has not be
            // included
            double variance = 0.0;
            for (int j = 0; j < nSpecdens; j++) {
                double covjk = 0.0;
                for (int i = 0; i < nReplicates; i++) {
                    double weight = scores[i].getWeights()[j];
                    double weightMean = weightMeans[j];
                    double parameter = parameters[k][i];
                    covjk += (weight - weightMean) * (parameter - mean);
                }
                covjk /= nReplicates;
                variance += Math.pow(covjk, 2.0);
            }
            double error = Math.sqrt(variance);
            orderPar = orderPar.set(name, mean, error);
        }
    }

    public double scoreBootstrap(Map<String, MolDataValues> molDataRes, MFModelIso model, double[] pars,
                                 BootstrapAggregator bootstrapAggregator, List<Integer> iRepList, int iStart, int nReplicates,
                                 boolean localFitTau, double localTauFraction) {
        double rssSum = 0.0;
        int startPar = localFitTau ? 0 : 1;
        double[] pars2 = new double[pars.length - startPar];

        System.arraycopy(pars, startPar, pars2, 0, pars2.length);
        for (int iRep = 0; iRep < nReplicates; iRep++) {
            int bootStrapSet = iRepList.get(iRep + iStart);
            for (var molData : molDataRes.values()) {
                molData.setBootstrapAggregator(bootstrapAggregator);
                molData.setBootstrapSet(bootStrapSet);
            }
            model.setTauFraction(localTauFraction);
            rssSum += scoreModel(molDataRes, pars2);
        }
        return rssSum / nReplicates;
    }

    Score tryModel(
        Map<String, MolDataValues> molDataRes,
        MFModelIso model,
        double localTauFraction,
        boolean localFitTau,
        Random random
    ) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambdaS(lambdaS);
        relaxFit.setLambdaTauF(lambdaTauF);
        relaxFit.setLambdaTauS(lambdaTauS);
        relaxFit.setUseLambda(useLambda);
        relaxFit.setFitJ(fitJ);

        model.setTauFraction(localTauFraction);

        double[] start = model.getStart();
        double[] lower = model.getLower();
        if (useLambda) {
            if (!(model instanceof MFModelIso2sf)) {
                throw new AssertionError(
                    "`useLambda` is set to `true` and the model type is not `MFModelIso2sf`. " +
                    "This should not be possible.");
            }
            // If regularization is used, all parameters should have lower bounds of 0.0
            lower = new double[lower.length];
        }
        double[] upper = model.getUpper();

        Optional<PointValuePair> result = relaxFit.fitResidueToModel(start, lower, upper);
        if (result.isPresent()) return relaxFit.score(result.get().getPoint(), true);
        else return null;
    }

    public void setTauFraction(double value) {
        this.tauFraction = value;
    }

    private void appendValueError(StringBuilder stringBuilder, Double val, Double err, String format) {
        if (val != null) {
            stringBuilder.append(String.format(format, val)).append("\t");
        } else {
            stringBuilder.append("\t\t");
        }
        if (err != null) {
            stringBuilder.append(String.format(format, err)).append("\t");
        } else {
            stringBuilder.append("\t\t");
        }
    }

    public void writeData(File file) throws IOException {
        try (FileWriter fileWriter = new FileWriter(file)) {
            var molDataValues = getData(false);
            String header = "residue\tfield\tr1\tr1_err\tr2\tr2_err\tnoe\tnoe_err\n";
            fileWriter.write(header);
            for (var molData : molDataValues.values()) {
                var data = molData.getData();
                for (var value : data) {
                    R1R2NOEDataValue dValue = (R1R2NOEDataValue) value;
                    StringBuilder stringBuilder = new StringBuilder();
                    stringBuilder.append(molData.specifier).append("\t");
                    stringBuilder.append(String.format("%.2f", dValue.relaxObj.getSF() / 1.0e6)).append("\t");
                    appendValueError(stringBuilder, dValue.R1, dValue.R1err, "%.3f");
                    appendValueError(stringBuilder, dValue.R2, dValue.R2err, "%.3f");
                    appendValueError(stringBuilder, dValue.NOE, dValue.NOEerr, "%.3f");
                    fileWriter.write(stringBuilder.toString());
                    fileWriter.write("\n");
                }
            }
        }
    }
}
