/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.rng.sampling.distribution.DirichletSampler;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;
import org.nmrfx.chemistry.*;
import org.nmrfx.chemistry.relax.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author brucejohnson
 */
public class FitR1R2NOEModel extends FitModel {
    Map<String, MolDataValues> molData = null;

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
            if (tau == null) {
                tau = estimateTau(molData).get("tau");
            }

            Map<String, ModelFitResult> results = testModels(molData, modelNames);
            return results;
        } else {
            throw new IllegalStateException("No relaxation data to analyze.  Need T1,T2 and NOE");
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
        if (tau == null) {
            tau = estimateTau(molData).get("tau");
        }
        AtomicInteger counts = new AtomicInteger();
        int n = molData.entrySet().size();
        Map<String, ModelFitResult> results = new HashMap<>();
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
                if (bootstrapMode != BootstrapMode.PARAMETRIC) {
                    Optional<ModelFitResult> result = testModelsWithBootstrapAggregation(orderParSetMap, resData, key, modelNames, random);
                    result.ifPresent(o -> results.put(key, o));
                } else {
                    Optional<ModelFitResult> result = testModels(orderParSetMap, resData, key, modelNames, random);
                    result.ifPresent(o -> results.put(key, o));
                }
            }
            int iCount = counts.incrementAndGet();
        });
        return results;
    }

    private OrderPar makeOrderPar(OrderParSet orderParSet, MolDataValues resData, Map<String, MolDataValues> molDataRes, String key,
                                  Score bestScore, MFModelIso bestModel, double[][] repData) {
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

    public Optional<ModelFitResult> testModels(Map<String, OrderParSet> orderParSetMap, MolDataValues resData, String key, List<String> modelNames, Random random) {
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        molDataRes.put(key, resData);

        MFModelIso bestModel = null;
        Score bestScore = null;
        double lowestAIC = Double.MAX_VALUE;
        boolean localFitTau;
        double localTauFraction;
        if (fitTau(molDataRes)) {
            localFitTau = true;
            localTauFraction = tauFraction;
        } else {
            localTauFraction = 0.0;
            localFitTau = false;
        }
        for (var modelName : modelNames) {
            MFModelIso model = MFModelIso.buildModel(modelName,
                    localFitTau, tau, localTauFraction, fitExchange);
            resData.setTestModel(model);
            Score score = tryModel(molDataRes, model, localTauFraction, localFitTau, random);
            OrderParSet orderParSet = orderParSetMap.get("order_parameter_list_"+ modelName);
            double[][] repData = null;
            if (nReplicates > 2) {
                double[] pars = score.getPars();
                repData = replicates(molDataRes, model, localTauFraction, localFitTau, pars, random);
                OrderPar orderPar = makeOrderPar(orderParSet, resData, molDataRes, key, score, model, repData);
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
            OrderPar orderPar = makeOrderPar(orderParSet, resData, molDataRes, key, bestScore, bestModel, replicateData);
            ModelFitResult modelFitResult = new ModelFitResult(orderPar, replicateData, null);
            result = Optional.of(modelFitResult);
        }
        return result;
    }

    public Optional<ModelFitResult> testModelsWithBootstrapAggregation(Map<String, OrderParSet> orderParSetMap, MolDataValues resData, String key, List<String> modelNames, Random random) {
        Optional<ModelFitResult> result = Optional.empty();
        Map<String, MolDataValues> molDataRes = new TreeMap<>();
        molDataRes.put(key, resData);

        boolean localFitTau;
        double localTauFraction;
        if (fitTau(molDataRes)) {
            localFitTau = true;
            localTauFraction = tauFraction;
        } else {
            localTauFraction = 0.0;
            localFitTau = false;
        }
        int maxPars = 5;
        int nExp = resData.dataValues.size();
        int nJ = nExp * 3;
        BootstrapAggregator bootstrapAggregator = new BootstrapAggregator(nExp);

        List<Integer> iRepList = null;
        DirichletSampler dirichlet = null;
        if (bootstrapMode == BootstrapMode.BAYESIAN) {
            dirichlet = DirichletSampler.symmetric(getRandomSource(), nJ, 4.0);
        } else {
            nReplicates = Math.min(nReplicates, bootstrapAggregator.getN());
            iRepList = IntStream.range(0, bootstrapAggregator.getN()).boxed().collect(Collectors.toList());
            Collections.shuffle(iRepList);
        }
        double[][] replicateData = new double[maxPars][nReplicates];
        MFModelIso[] bestModels = new MFModelIso[nReplicates];
        Score[] bestScores = new Score[nReplicates];

        int[] totalCounts = new int[nJ];
        for (int iRep = 0; iRep < nReplicates; iRep++) {
            if (cancelled.get()) {
                return result;
            }
            if (bootstrapMode == BootstrapMode.BAYESIAN) {
                double[] weights = dirichlet.sample();
                scaleWeights(weights);
                for (var molData : molDataRes.values()) {
                    molData.weight(weights);
                    molData.clearBootStrapSet();
                }
            } else {
                int bootStrapSet = iRepList.get(iRep);
                for (var molData : molDataRes.values()) {
                    molData.setBootstrapAggregator(bootstrapAggregator);
                    molData.setBootstrapSet(bootStrapSet);
                    molData.weight(null);
                }
                BootstrapAggregator.incrCounts(totalCounts, bootstrapAggregator.getY(bootStrapSet));
            }
            List<Score> scores = new ArrayList<>();
            List<MFModelIso> models = new ArrayList<>();
            for (var modelName : modelNames) {
                MFModelIso model = MFModelIso.buildModel(modelName,
                        localFitTau, tau, localTauFraction, fitExchange);
                resData.setTestModel(model);
                Score score = tryModel(molDataRes, model, localTauFraction, localFitTau, random);
                scores.add(score);
                models.add(model);
            }
            int iBest = bestScore(scores);
            if (iBest != -1) {
                Score bestScore = scores.get(iBest);
                MFModelIso bestModel = models.get(iBest);
                double[] repPars = bestScore.getPars();
                bestModels[iRep] = bestModel;
                bestScores[iRep] = bestScore;
                double[] pars = bestModel.getStandardPars(repPars);

                for (int iPar = 0; iPar < pars.length; iPar++) {
                    replicateData[iPar][iRep] = pars[iPar];
                }
            }
        }
        for (int i=0;i< totalCounts.length;i++) {
            totalCounts[i] /= nReplicates;
        }
        for (var molData : molDataRes.values()) {
            molData.clearBootStrapSet();;
        }
        MFModelIso bestModel = MFModelIso.buildModel("2sf",
                localFitTau, tau, localTauFraction, fitExchange);

        if (bestModels[0] != null) {
            String[] parNames = {"Tau_e", "Sf2", "Tau_f", "Ss2", "Tau_s"};

            ResonanceSource resSource = new ResonanceSource(resData.atom);
            double rssSum = 0.0;
            for (int i = 0; i < nReplicates; i++) {
                rssSum += bestScores[i].rss;
            }
            double rss = rssSum /= nReplicates;
            OrderParSet orderParSet = orderParSetMap.get("order_parameter_list_best");
            OrderPar orderPar = new OrderPar(orderParSet, resSource, rss, bestScores[0].nValues, parNames.length, bestModel.getName());
            double[][] cov = new double[nJ][parNames.length];
            double[] bestPars = new double[parNames.length];
            for (int iPar = 0; iPar < parNames.length; iPar++) {
                String parName = parNames[iPar];
                Double parError = null;
                DescriptiveStatistics sumStat = new DescriptiveStatistics(replicateData[iPar]);
                double parValue = useMedian ? sumStat.getPercentile(50.0) : sumStat.getMean();
                bestPars[iPar] = parValue;
                if (bootstrapMode == BootstrapMode.BAYESIAN) {
                    parError = sumStat.getStandardDeviation();
                } else {
                    for (int j = 0; j < nJ; j++) {
                        for (int i = 0; i < nReplicates; i++) {
                            int[] y = bootstrapAggregator.getY(iRepList.get(i));
                            cov[j][iPar] += (y[j] - totalCounts[j]) * (replicateData[iPar][i] - parValue);
                        }
                        cov[j][iPar] /= nReplicates;
                    }

                    double sum = 0.0;
                    for (int j = 0; j < nJ; j++) {
                        sum += cov[j][iPar] * cov[j][iPar];
                    }
                    parError = Math.sqrt(sum / nJ);
                }
                orderPar = orderPar.set(parName, parValue, parError);
            }
            orderPar = orderPar.set("model", (double) bestModel.getNumber(), null);
            Atom atom = resData.atom;
            atom.addOrderPar(orderParSet, orderPar);
            double[][] jValues = resData.getJValues();
            SpectralDensity spectralDensity = new SpectralDensity(key, jValues);
            atom.addSpectralDensity(key, spectralDensity);
            resData.setTestModel(bestModel);
            Double validationScore = null;
            if (calcValidation) {
                if (bootstrapMode == BootstrapMode.BAYESIAN) {
                    validationScore = scoreBayesian(molDataRes, bestModel, bestPars, dirichlet,
                            nReplicates, nReplicates, localFitTau, localTauFraction);
                } else {
                    if ((nReplicates + nReplicates) <= iRepList.size()) {
                        validationScore = scoreBootstrap(molDataRes, bestModel, bestPars, bootstrapAggregator, iRepList,
                                nReplicates, nReplicates, localFitTau, localTauFraction);
                    }
                }
            }
            ModelFitResult modelFitResult = new ModelFitResult(orderPar, replicateData, validationScore);
            result = Optional.of(modelFitResult);
        }
        return result;
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

    Score tryModel(Map<String, MolDataValues> molDataRes, MFModelIso model, double localTauFraction, boolean localFitTau, Random random) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambdaS(lambdaS);
        relaxFit.setLambdaTau(lambdaTau);
        relaxFit.setUseLambda(useLambda);
        relaxFit.setFitJ(fitJ);
        model.setTauFraction(localTauFraction);
        double[] start = model.getStart();
        double[] lower = model.getLower();
        double[] upper = model.getUpper();
        double[] keepStart = start.clone();
        int nTries = 3;
        PointValuePair best = null;
        for (int i = 0; i < nTries; i++) {
            Optional<PointValuePair> fitResultOpt = relaxFit.fitResidueToModel(start, lower, upper);
            if (fitResultOpt.isPresent() && ((i == 0) || (fitResultOpt.get().getValue() < best.getValue()))) {
                best = fitResultOpt.get();
            }
            for (int j = 0; j < start.length; j++) {
                start[j] = keepStart[j] + random.nextGaussian() * 0.1 * (upper[j] - lower[j]);
            }
        }
        if (best != null) {
            var score = relaxFit.score(best.getPoint(), true);
            return score;
        } else {
            return null;
        }
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
