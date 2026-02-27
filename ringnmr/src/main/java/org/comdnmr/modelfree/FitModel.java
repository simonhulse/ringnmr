package org.comdnmr.modelfree;

import javafx.beans.property.ReadOnlyObjectProperty;
import javafx.concurrent.Service;
import javafx.concurrent.Task;
import javafx.concurrent.Worker;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.DirichletSampler;
import org.apache.commons.rng.simple.RandomSource;
import org.comdnmr.BasicFitter;
import org.comdnmr.eqnfit.ParValueInterface;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.util.ProcessingStatus;

import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.Function;

public abstract class FitModel implements BasicFitter {
    public static UniformRandomProvider rng = null;

    Double tau;
    boolean fitTau = false;
    boolean fitJ = false;
    BootstrapMode bootstrapMode = BootstrapMode.PARAMETRIC;
    boolean bayesian = false;
    boolean fitExchange = false;
    double tauFraction = 0.25;
    double lambdaS = 0.0;
    double lambdaTau = 0.0;
    boolean useLambda = false;
    double t2Limit = 0.0;
    int nReplicates = 0;
    private final FitResidues fitResidues = new FitResidues();
    final ReadOnlyObjectProperty<Worker.State> stateProperty = fitResidues.worker.stateProperty();
    private boolean isProcessing = false;
    List<String> modelNames = new ArrayList<>();
    String searchKey = null;
    AtomicBoolean cancelled = new AtomicBoolean(false);
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

    @Override
    public List<ParValueInterface> guessPars(String eqn) {
        return List.of();
    }

    @Override
    public double[] getSimXDefaults() {
        return new double[0];
    }


    public static UniformRandomProvider getRandomSource() {
        if (rng == null) {
            rng = RandomSource.XO_RO_SHI_RO_128_PP.create();
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

    public double scoreBayesian(Map<String, MolDataValues> molDataRes, MFModelIso model, double[] pars,
                                DirichletSampler dirichlet, int iStart, int nReplicates,
                                boolean localFitTau, double localTauFraction) {
        double rssSum = 0.0;
        int startPar = localFitTau ? 0 : 1;
        double[] pars2 = new double[pars.length - startPar];

        System.arraycopy(pars, startPar, pars2, 0, pars2.length);
        for (int iRep = 0; iRep < nReplicates; iRep++) {
            double[] weights = dirichlet.sample();
            scaleWeights(weights);
            for (var molData : molDataRes.values()) {
                molData.weight(weights);
            }
            model.setTauFraction(localTauFraction);
            rssSum += scoreModel(molDataRes, pars2);
        }
        return rssSum / nReplicates;
    }

     double scoreModel(Map<String, MolDataValues> molDataRes, double[] pars) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambdaS(lambdaS);
        relaxFit.setLambdaTau(lambdaTau);
        relaxFit.setUseLambda(useLambda);
        relaxFit.setFitJ(fitJ);
        var score = relaxFit.score(pars, true);
        return score.rss;
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

    public abstract Map<String, ModelFitResult> testIsoModel();

    double[][] replicates(Map<String, MolDataValues> molDataRes,
                          MFModelIso bestModel, double localTauFraction,
                          boolean localFitTau, double[] pars, Random random) {
        double[][] repData = new double[pars.length][nReplicates];
        for (int iRep = 0; iRep < nReplicates; iRep++) {
            Score score2 = fitReplicate(molDataRes, bestModel, localTauFraction, localFitTau, pars, random);
            double[] repPars = score2.getPars();
            for (int iPar = 0; iPar < pars.length; iPar++) {
                repData[iPar][iRep] = repPars[iPar];
            }
        }
        return repData;
    }

    Score fitReplicate(Map<String, MolDataValues> molDataRes, MFModelIso model,
                       double localTauFraction, boolean localFitTau, double[] pars, Random random) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(molDataRes);
        relaxFit.setLambdaS(lambdaS);
        relaxFit.setLambdaTau(lambdaTau);
        relaxFit.setUseLambda(useLambda);
        relaxFit.setFitJ(fitJ);
        Map<String, MolDataValues> molDataMap = relaxFit.genBootstrap(random, model, pars);
        relaxFit.setRelaxData(molDataMap);

        model.setTauFraction(localTauFraction);
        double[] lower = model.getLower();
        double[] upper = model.getUpper();
        Optional<PointValuePair> fitResultOpt = relaxFit.fitResidueToModel(pars, lower, upper);
        if (fitResultOpt.isPresent()) {
            return relaxFit.score(fitResultOpt.get().getPoint(), true);
        }
        return null;
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

    public void setLambdaS(double value) {
        this.lambdaS = value;
    }

    public void setLambdaTau(double value) {
        this.lambdaTau = value;
    }

    public void setUseLambda(boolean value) {
        this.useLambda = value;
    }

    public boolean useLambda() {
        return useLambda && (lambdaS > 1.0e-8 || lambdaTau > 1.0e-8);
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

    public boolean fitTau(Map<String, MolDataValues> molDataRes) {
        boolean localFitTau;

        if (overT2Limit(molDataRes, t2Limit)) {
            localFitTau = fitTau;
        } else {
            localFitTau = false;
        }
        return localFitTau;
    }

    boolean overT2Limit(Map<String, MolDataValues> molDataRes, double limit) {
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
        updateStatus("Done");
        try {
        } catch (IllegalArgumentException iAE) {
        }
    }

    public void fitAll(Task task) {
        testIsoModel();
        int nFit = 0;
        int nGroups = 1;
        nFit++;
        if (task != null) {
            updateProgress((1.0 * nFit) / nGroups);
        }
    }
}

