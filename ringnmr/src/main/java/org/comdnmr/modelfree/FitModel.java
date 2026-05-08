package org.comdnmr.modelfree;

import javafx.beans.property.ReadOnlyObjectProperty;
import javafx.concurrent.Service;
import javafx.concurrent.Task;
import javafx.concurrent.Worker;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.comdnmr.BasicFitter;
import org.comdnmr.eqnfit.ParValueInterface;
import org.comdnmr.util.ProcessingStatus;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;

import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.relax.OrderParSet;

public abstract class FitModel implements BasicFitter {
    public static UniformRandomProvider rng = null;

    protected FitSpec fitSpec = new ConventionalFitSpec.Builder()
        .modelNames(List.of("1", "1f", "1s", "2s", "2sf"))
        .build();

    protected StructureValues molData = null;

    // >>>>>>>>>>>>>>>>>>>>>>>>
    // TODO: transition to not needing these attributes. It is all contained
    // within FitSpec
    Double tau;
    boolean fitTau = false;
    boolean fitJ = false;
    BootstrapMode bootstrapMode = BootstrapMode.PARAMETRIC;
    boolean bayesian = false;
    boolean fitExchange = false;
    double tauFraction = 0.25;
    double lambdaS2F = 0.0;
    double lambdaS2S = 0.0;
    double lambdaTauF = 0.0;
    double lambdaTauS = 0.0;
    boolean useLambda = false;
    double t2Limit = 0.0;
    int nReplicates = 0;
    // <<<<<<<<<<<<<<<<<<<<<<<<

    private final FitResidues fitResidues = new FitResidues();
    final ReadOnlyObjectProperty<Worker.State> stateProperty = fitResidues.worker.stateProperty();
    private boolean isProcessing = false;
    List<String> modelNames = new ArrayList<>();
    String searchKey = null;
    AtomicBoolean cancelled = new AtomicBoolean(false);
    int lastFitCount = 0;
    boolean useMedian = false;
    boolean calcValidation = false;

    Function<Double, Double> updaterFunction;
    Function<ProcessingStatus, Double> statusFunction;

    public enum BootstrapMode {
        PARAMETRIC,
        AGGREGATE,
        BAYESIAN
    }

    @Override
    public double rms(double[] pars) {
        return 0;
    }

    @Override
    public void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues) {

    }

    public void setFitSpec(FitSpec fitSpec) { this.fitSpec = fitSpec; }

    public void setData(StructureValues molData) { this.molData = molData; }

    /**
     * Loads relaxation data from the active molecule. Called when {@link #molData}
     * is null or empty at the start of {@link #testIsoModel()}.
     */
    protected abstract StructureValues loadData();

    /**
     * Sets the overall tumbling correlation time on {@link #fitSpec} when it has
     * not yet been computed, using the estimate provided by {@link #molData}.
     */
    protected void setTauMFromData() {
        fitSpec.setTauM(molData.estimateTau().get("tau"));
    }

    /**
     * Returns the error message thrown when no relaxation data are available.
     */
    protected abstract String getNoDataMessage();

    @Override
    public List<ParValueInterface> guessPars(String eqn) {
        return List.of();
    }

    @Override
    public double[] getSimXDefaults() {
        return new double[0];
    }


    public static UniformRandomProvider getRandomSource() {
        return getRandomSource(false);
    }

    public static UniformRandomProvider getRandomSource(boolean init) {
        if (init || rng == null) {
            final int[] seed = new int[] { 196, 9, 0, 226 };
            rng = RandomSource.XO_RO_SHI_RO_128_PP.create(seed);
        }
        return rng;
    }


    public void setup(String searchKey, List<String> modelNames) {
        this.searchKey = searchKey;
        this.modelNames.clear();
        this.modelNames.addAll(modelNames);
    }

    void scaleWeights(double[] weights) {
        for (int i = 0; i < weights.length; i++) {
            weights[i] = weights[i] * weights.length;
        }

    }

    int bestScore(List<Score> scores) {
        double lowestAIC = Double.MAX_VALUE;
        int iBest = -1;
        int i = 0;
        for (var score : scores) {
            if (score.aic() < lowestAIC) {
                lowestAIC = score.aic();
                iBest = i;
            }
            i++;
        }
        return iBest;
    }

    public void updaters(Function<Double, Double> updaterFunction, Function<ProcessingStatus, Double> statusFunction) {
        this.updaterFunction = updaterFunction;
        this.statusFunction = statusFunction;
    }

    public Map<String, ModelFitResult> testIsoModel() {
        if ((molData == null) || (molData.isEmpty())) {
            molData = loadData();
        }
        if ((searchKey != null) && molData.containsKey(searchKey)) {
            var keepVal = molData.get(searchKey);
            molData.clear();
            molData.put(searchKey, keepVal);
        }
        if (!molData.isEmpty()) {
            if (fitSpec.tauMNeedsComputing()) setTauMFromData();
            AtomicInteger counts = new AtomicInteger();
            int n = molData.size();
            MoleculeBase moleculeBase = MoleculeFactory.getActive();
            Map<String, OrderParSet> orderParSetMap = new ConcurrentHashMap<>();
            orderParSetMap.putAll(moleculeBase.orderParSetMap());
            Map<String, ModelFitResult> results = new ConcurrentHashMap<>();
            new ArrayList<>(molData.entrySet())
                .parallelStream()
                .forEach(residue -> {
                    updateProgress((double) counts.get() / n);
                    if (cancelled.get()) return;
                    String key = residue.getKey();
                    MolDataValues<? extends RelaxDataValue> data = residue.getValue();
                    if (!data.getData().isEmpty()) results.put(key, fitSpec.fit(key, data, orderParSetMap));
                    counts.incrementAndGet();
                });
            moleculeBase.orderParSetMap().putAll(orderParSetMap);
            return results;
        } else {
            throw new IllegalStateException(getNoDataMessage());
        }
    }

    public void setUseMedian(boolean value) {
        useMedian = value;
    }

    public void setCalcValidation(boolean value) {
        calcValidation = value;
    }

    public Double getTau() {
        return tau;
    }

    public void setTau(Double value) {
        tau = value;
    }

    public void setFitTau(boolean value) {
        fitTau = value;
    }

    public void setLambdaS2F(double value) {
        this.lambdaS2F = value;
    }

    public void setLambdaS2S(double value) {
        this.lambdaS2S = value;
    }

    public void setLambdaTauF(double value) {
        this.lambdaTauF = value;
    }

    public void setLambdaTauS(double value) {
        this.lambdaTauS = value;
    }

    public void setUseLambda(boolean value) {
        this.useLambda = value;
    }

    public boolean useLambda() {
        return useLambda && (lambdaS2F > 1.0e-8 || lambdaS2S > 1.0e-8 || lambdaTauF > 1.0e-8 || lambdaTauS > 1.0e-8);
    }

    public void setFitJ(boolean value) {
        this.fitJ = value;
    }

    public void setBootstrapMode(BootstrapMode value) {
        this.bootstrapMode = value;
    }

    public void setNReplicates(int value) {
        this.nReplicates = value;
    }

    public int getNReplicates() {
        return nReplicates;
    }

    public void setT2Limit(double value) {
        this.t2Limit = value;
    }

    public void setTauFraction(double value) {
        tauFraction = value;
    }

    public boolean fitTau(Map<String, MolDataValues<? extends RelaxDataValue>> molDataRes) {
        boolean localFitTau;

        if (overT2Limit(molDataRes, t2Limit)) {
            localFitTau = fitTau;
        } else {
            localFitTau = false;
        }
        return localFitTau;
    }

    boolean overT2Limit(Map<String, MolDataValues<? extends RelaxDataValue>> molDataRes, double limit) {
        return limit < 1.0e-6 || molDataRes.values().stream().anyMatch(v -> v.getData().stream().anyMatch(d -> d.R2 > limit));
    }

    private class FitResidues {

        String script;
        public Worker<Integer> worker;

        private FitResidues() {
            worker = new Service<Integer>() {

                @Override
                protected Task createTask() {
                    return new Task() {
                        @Override
                        protected Object call() {
                            cancelled.set(false);
                            updateStatus("Start processing");
                            updateTitle("Start Processing");
                            fitAll(this);
                            return 0;
                        }
                    };
                }
            };

            ((Service<Integer>) worker).setOnSucceeded(event -> {
                finishProcessing();
                setProcessingOff();
            });
            ((Service<Integer>) worker).setOnCancelled(event -> {
                cancelled.set(true);
                setProcessingOff();
                setProcessingStatus("cancelled", false);
            });
            ((Service<Integer>) worker).setOnFailed(event -> {
                setProcessingOff();
                final Throwable exception = worker.getException();
                setProcessingStatus(exception.getMessage(), false, exception);

            });

        }
    }


    public void fitResidues() {
        setProcessingOn();
        updateProgress(0.0);
        ((Service) fitResidues.worker).restart();
    }

    public void updateProgress(Double f) {
        if (updaterFunction != null) {
            updaterFunction.apply(f);
        }
    }

    public void updateStatus(String s) {
        setProcessingStatus(s, true);
    }

    public void setProcessingStatus(String s, boolean ok) {
        setProcessingStatus(s, ok, null);
    }

    public void setProcessingStatus(String s, boolean ok, Throwable throwable) {
        if (statusFunction != null) {
            ProcessingStatus status = new ProcessingStatus(s, ok, throwable);
            statusFunction.apply(status);
        }
    }

    public void clearProcessingTextLabel() {
        ProcessingStatus status = new ProcessingStatus("");
        statusFunction.apply(status);
    }

    public void haltFit() {
        fitResidues.worker.cancel();
    }

    synchronized void setProcessingOn() {
        isProcessing = true;
    }

    synchronized void setProcessingOff() {
        isProcessing = false;
    }

    boolean isProcessing() {
        return isProcessing;
    }

    void finishProcessing() {
        updateStatus("Done: fit " + lastFitCount + " residues");
    }

    public void fitAll(Task task) {
        lastFitCount = testIsoModel().size();
        if (task != null) {
            updateProgress(1.0);
        }
    }
}

