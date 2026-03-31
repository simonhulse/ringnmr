package org.comdnmr.modelfree;

import java.util.*;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.SpectralDensity;


// TODO: if tauM is not specified, it must be estimated using the data.
// Including the Map<String, MolDataValues> into the build method will allow
// this to be computed without the risk if an IllegalStateException
// (see getTauM below). estimateTau() is currently defined in FitR1R2NOEModel. Makes more sense to me to make this a method of this class
public abstract class FitSpec {

    protected double tauM;
    protected boolean tauMNeedsComputing;
    protected final boolean fitTauM;
    protected final boolean fitJ; // TODO: Currently not supported
    protected final BootstrapMode bootstrapMode;
    protected final boolean fitExchange; // TODO: Currently not supported
    protected final double tauFraction;
    protected final double t2Limit;
    protected final int nReplicates;

    private static final Map<String, Class<? extends FitSpec>> CLASSES = new LinkedHashMap<>();
    static {
        CLASSES.put("Conventional", ConventionalFitSpec.class);
        CLASSES.put("Bootstrap Aggregation", BaggingFitSpec.class);
        CLASSES.put("Regularization", RegularizationFitSpec.class);
    }

    protected FitSpec(Builder<?> builder) {
        if (builder.tauM.isEmpty()) {
            this.tauM = 0.0;
            tauMNeedsComputing = true;
        } else {
            this.tauM = builder.tauM.get();
            tauMNeedsComputing = false;
        }
        this.fitTauM = builder.fitTauM;
        this.fitJ = builder.fitJ;
        this.bootstrapMode = builder.bootstrapMode;
        this.fitExchange = builder.fitExchange;
        this.tauFraction = builder.tauFraction;
        this.t2Limit = builder.t2Limit;
        this.nReplicates = builder.nReplicates;
    }

    public abstract ModelFitResult fit(String key, MolDataValues data);

    public static List<String> getNames() { return new ArrayList<>(CLASSES.keySet()); }

    public static Class<? extends FitSpec> getFitSpec(String methodName) {
        if (CLASSES.containsKey(methodName)) {
            return CLASSES.get(methodName);
        } else {
            throw new IllegalArgumentException(
                String.format(
                    "`methodName` should be on of: %s",
                    String.join(", ", getNames())
                )
            );
        }
    }

    protected static Score runFit(RelaxFit relaxFit, MFModelIso model) {
        double[] start = model.getStart();
        double[] lower = model.getLower();
        double[] upper = model.getUpper();

        Optional<PointValuePair> result = relaxFit.fitResidueToModel(start, lower, upper);
        if (result.isEmpty()) {
            throw new RuntimeException("Could not generate fit result.");
        }
        return relaxFit.score(result.get().getPoint(), true);
    }

    public boolean tauMNeedsComputing() {
        return tauMNeedsComputing;
    }

    public void setTauM(double tauM) {
        this.tauM = tauM;
        tauMNeedsComputing = false;
    }

    public double getTauM() {
        if (tauMNeedsComputing) {
            throw new IllegalStateException("tauM has not be set!");
        }
        return tauM;
    }

    public boolean fitTauM(MolDataValues data) {
        return fitTauM && data.getData().stream().anyMatch(value -> value.R2 > t2Limit);
    }

    public enum BootstrapMode {
        PARAMETRIC,
        NONPARAMETRIC,
        BAYESIAN
    }

    public BootstrapSampler getBootstrapSampler(MolDataValues data) {
        return switch (bootstrapMode) {
            case PARAMETRIC -> new ParametricSampler(data);
            case NONPARAMETRIC -> new NonparametricSampler(data);
            case BAYESIAN -> new BayesianSampler(data);
        };
    }

    protected RelaxFit initRelaxFit(String key, MolDataValues data) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setFitJ(fitJ);
        relaxFit.setRelaxData(key, data);
        return relaxFit;
    }

    private double[] computeMeans(double[][] array) {
        int nSamples = array.length;
        int nDataPoints = array[0].length;

        double[] means = new double[nDataPoints];
        for (int j = 0; j < nDataPoints; j++) {
            double sum = 0.0;
            for (int i = 0; i < nSamples; i++) {
                sum += array[i][j];
            }
            double mean = sum / nSamples;
            means[j] = mean;
        }
        return means;
    }

    private double[] computeStdevs(double[][] array, Optional<double[]> meansOpt) {
        int nSamples = array.length;
        int nDataPoints = array[0].length;

        double[] means;
        if (meansOpt.isEmpty()) means = computeMeans(array);
        else means = meansOpt.get();
        double[] stdevs = new double[nDataPoints];
        for (int j = 0; j < nDataPoints; j++) {
            double sum = 0.0;
            double mean = means[j];
            for (int i = 0; i < nSamples; i++) {
                sum += Math.pow(mean - array[i][j], 2.0);
            }
            double stdev = Math.sqrt(sum / (nSamples - 1));
            stdevs[j] = stdev;
        }
        return stdevs;
    }

    private Pair<double[], double[]> computeStatistics(double[][] parameters, double[][] weights) {
        return switch (bootstrapMode) {
            case PARAMETRIC -> computeStatisticsParametric(parameters);
            case NONPARAMETRIC, BAYESIAN -> computeStatisticsBootstrapping(parameters, weights);
        };
    }

    private Pair<double[], double[]> computeStatisticsParametric(double[][] parameters) {
        double[] parameterMeans = computeMeans(parameters);
        double[] parameterErrors = computeStdevs(parameters, Optional.of(parameterMeans));
        return Pair.of(parameterMeans, parameterErrors);
    }

    private Pair<double[], double[]> computeStatisticsBootstrapping(double[][] parameters, double[][] weights) {
        int nParameters = parameters[0].length;
        int nDataPoints = weights[0].length;

        // Compute weight means. Should be 1.0 for every datapoint for
        // Nonparametric, and ~1.0 for Bayesian
        double[] parameterMeans = computeMeans(parameters);
        double[] weightMeans = computeMeans(weights);

        double[] parameterErrors = new double[nParameters];
        for (int k = 0; k < nParameters; k++) {
            double parameterMean = parameterMeans[k];
            double variance = 0.0;
            for (int j = 0; j < nDataPoints; j++) {
                double covjk = 0.0;
                double weightMean = weightMeans[j];
                for (int i = 0; i < nReplicates; i++) {
                    double weight = weights[i][j];
                    double parameter = parameters[i][k];
                    covjk += (weight - weightMean) * (parameter - parameterMean);
                }
                covjk /= nReplicates;
                variance += Math.pow(covjk, 2.0);
            }
            double parameterError = Math.sqrt(variance);
            parameterErrors[k] = parameterError;
        }
        return Pair.of(parameterMeans, parameterErrors);
    }


    protected OrderPar makeOrderParSet(
        OrderParSet orderParSet,
        MolDataValues data,
        String key,
        Score score,
        MFModelIso model,
        double[][] parameters,
        double[][] weights
    ) {
        Pair<double[], double[]> parameterStats = computeStatistics(parameters, weights);
        double[] parameterMeans = parameterStats.getLeft();
        double[] parameterErrors = parameterStats.getRight();
        int nParameters = parameterMeans.length;

        ResonanceSource resSource = new ResonanceSource(data.atom);
        Atom atom = data.atom;
        List<String> parameterNames = model.getParNames();
        OrderPar orderPar = new OrderPar(orderParSet, resSource, score.rss, score.nValues, score.nPars, model.getName());

        for (int k = 0; k < nParameters; k++) {
            String parameterName = parameterNames.get(k);
            double parameterMean = parameterMeans[k];
            double parameterError = parameterErrors[k];
            orderPar = orderPar.set(parameterName, parameterMean, parameterError);
        }

        // If tauM is fixed, set it to the fixed value with zero error
        if (!model.fitTau())  orderPar = orderPar.set("Tau_e", model.getTau(), 0.0);

        orderPar = orderPar.set("model", (double) model.getNumber(), null);
        atom.addOrderPar(orderParSet, orderPar);

        SpectralDensity spectralDensity = new SpectralDensity(key, data.getJValues());
        atom.addSpectralDensity(key, spectralDensity);

        return orderPar;
    }

    public static abstract class Builder<T extends Builder<T>> {
        private Optional<Double> tauM = Optional.empty();
        private boolean fitTauM = false;
        private boolean fitJ = true;
        private BootstrapMode bootstrapMode = BootstrapMode.PARAMETRIC;
        private boolean fitExchange = false; // Currently not supported
        private double tauFraction = 0.25;
        private double t2Limit = 0.0;
        // 25 is smaller than 27, the maximum number of iterates possible for a
        // 2-field dataset when using nonparametric bootstrapping.
        // It therefore will not lead to an unexpected error if the user simply
        // uses the default setup
        private int nReplicates = 25;

        @SuppressWarnings("unchecked")
        private T self() {
            return (T) this;
        }

        public T tauM(Optional<Double> tauM) {
            this.tauM = tauM;
            return self();
        }

        public T fitTauM(boolean fitTauM) {
            this.fitTauM = fitTauM;
            return self();
        }

        // TODO: implement
        public T fitJ(boolean fitJ) {
            // this.fitJ = fitJ;
            this.fitJ = true;
            return self();
        }

        public T bootstrapMode(BootstrapMode bootstrapMode) {
            this.bootstrapMode = bootstrapMode;
            return self();
        }

        // TODO: implement
        public T fitExchange(boolean fitExchange) {
            // this.fitExchange = fitExchange;
            this.fitExchange = false;
            return self();
        }

        public T tauFraction(double tauFraction) {
            this.tauFraction = tauFraction;
            return self();
        }

        public T t2Limit(double t2Limit) {
            this.t2Limit = t2Limit;
            return self();
        }

        public T nReplicates(int nReplicates) {
            this.nReplicates = nReplicates;
            return self();
        }

        public abstract FitSpec build();
    }
}
