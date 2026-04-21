package org.comdnmr.modelfree;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.comdnmr.modelfree.models.MFModelIso;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.SpectralDensity;

/**
 * Abstract base class that defines the specification and execution strategy for
 * model-free analysis fitting of NMR relaxation data.
 *
 * <p>Concrete subclasses (e.g. {@link ConventionalFitSpec}, {@link BaggingFitSpec},
 * {@link RegularizationFitSpec}) implement specific fitting methodologies while
 * sharing common configuration such as the overall tumbling correlation time
 * (tauM), bootstrap resampling mode, TOML serialization, etc.</p>
 *
 * <p>Instances are constructed via the nested {@link Builder} using a
 * curiously-recurring template pattern (CRTP) to support subclass-specific
 * builder methods while preserving fluent chaining.</p>
 *
 * <p><strong>Note on tauM estimation:</strong> If tauM is not specified at
 * construction time, it must be estimated from the data before fitting. The
 * {@link #tauMNeedsComputing()} and {@link #setTauM(double)} methods manage
 * this lifecycle. See also {@code FitR1R2NOEModel.estimateTau()} which
 * currently handles the estimation.</p>
 *
 * @see ConventionalFitSpec
 * @see BaggingFitSpec
 * @see RegularizationFitSpec
 */
public abstract class FitSpec {

    /** The type of moiety being analyzed (e.g. amide backbone or deuterated methyl). */
    protected MoietyType moietyType;

    /**
     * Overall tumbling correlation time in seconds. Zero when not yet estimated.
     * Access through {@link #getTauM()} which enforces the computation guard.
     */
    private double tauM;

    /**
     * Flag indicating whether {@link #tauM} still needs to be computed from data.
     * Managed internally by {@link #setTauM(double)} and the constructor.
     */
    private boolean tauMNeedsComputing;

    /** Whether tau_M should be included as a free parameter in the fit. */
    protected final boolean fitTauM;

    /** Whether spectral density J values should be fit directly. Not yet supported. */
    protected final boolean fitJ;

    /** The bootstrap resampling strategy used for error estimation. */
    protected final BootstrapMode bootstrapMode;

    /** Whether chemical-exchange contributions (Rex) should be fit. Not yet supported. */
    protected final boolean fitExchange;

    /**
     * Fractional bound on tau_M during fitting. If the initial tau_M estimate
     * is τ₀, the optimizer constrains tau_M to [τ₀·(1 − tauMFraction), τ₀·(1 + tauMFraction)].
     */
    protected final double tauMFraction;

    /**
     * R2 threshold: tau_M fitting is only attempted for residues with at least
     * one R2 value which exceeds this limit.
     */
    protected final double r2Limit;

    /** Number of bootstrap replicates to generate for error estimation. */
    protected final int nReplicates;

    /**
     * If {@code true}, use the median rather than the mean as the central
     * estimator when aggregating bootstrap replicates.
     *
     * <p><strong>Note:</strong> This setting is only meaningful for Bagging and
     * Nonparametric modes. The Conventional method uses the fit to the original
     * data for parameter estimates.</p>
     */
    protected final boolean useMedian;

    /**
     * Registry mapping human-readable method names to their corresponding
     * {@code FitSpec} subclass. Insertion order is preserved via
     * {@link LinkedHashMap}.
     */
    private static final Map<String, Class<? extends FitSpec>> CLASSES = new LinkedHashMap<>();
    static {
        CLASSES.put("Conventional", ConventionalFitSpec.class);
        CLASSES.put("Bootstrap Aggregation", BaggingFitSpec.class);
        CLASSES.put("Regularization", RegularizationFitSpec.class);
    }

    /**
     * Constructs a new {@code FitSpec} from the given builder.
     *
     * @param builder the builder containing all configuration values; must not be {@code null}
     */
    protected FitSpec(Builder<?> builder) {
        this.moietyType = builder.moietyType;
        if (builder.tauM.isEmpty()) {
            this.tauM = 0.0;
            this.tauMNeedsComputing = true;
        } else {
            this.tauM = builder.tauM.get();
            this.tauMNeedsComputing = false;
        }
        this.fitTauM = builder.fitTauM;
        this.fitJ = builder.fitJ;
        this.bootstrapMode = builder.bootstrapMode;
        this.fitExchange = builder.fitExchange;
        this.tauMFraction = builder.tauMFraction;
        this.r2Limit = builder.r2Limit;
        this.nReplicates = builder.nReplicates;
        this.useMedian = builder.useMedian;
    }

    /**
     * Builds a single {@link MFModelIso} by name, configured with the
     * current tau_M, tau_M fraction bounds, and exchange-fitting flag.
     *
     * @param name the model name
     * @param data relaxation data for the residue, used to determine whether
     *             tau_M should be fit for this residue
     * @return the configured model instance
     * @throws IllegalStateException if tau_M has not been set
     */
    protected MFModelIso getModel(String name, MolDataValues<? extends RelaxDataValue> data) {
        boolean fitTauM = fitTauM(data);
        String fullName = moietyType.getModelName(name);
        return MFModelIso.buildModel(fullName, fitTauM, getTauM(), tauMFraction, fitExchange);
    }


    // ── Abstract methods ────────────────────────────────────────────────

    /**
     * Performs the full model-free fit for a single residue.
     *
     * @param key            residue identifier
     * @param data           measured relaxation values for the residue
     * @param orderParSetMap mutable map to which the resulting {@link OrderParSet}
     *                       entries will be added
     * @return the result of the fitting process
     */
    public abstract ModelFitResult fit(String key, MolDataValues<?> data, Map<String, OrderParSet> orderParSetMap);

    /**
     * Returns the lower bounds for the model parameters. Subclasses may apply
     * strategy-specific adjustments (e.g. for regularization, all parameters
     * are bounded from below by 0).
     *
     * @param model the isotropic model-free model
     * @return array of lower-bound values, one per local parameter
     */
    protected double[] getLower(MFModelIso model) { return model.getLower(); }

    /**
     * Serializes this specification to a TOML-formatted string suitable for
     * writing to a configuration file.
     *
     * @return TOML representation of this fit specification
     */
    public abstract String toToml();

    // ── Static registry methods ─────────────────────────────────────────

    /**
     * Returns the list of registered fitting-method display names.
     *
     * @return an unmodifiable list of the registered method names
     */
    public static List<String> getNames() {
        return List.copyOf(CLASSES.keySet());
    }

    /**
     * Looks up the {@code FitSpec} subclass registered under the given name.
     *
     * @param methodName one of the names returned by {@link #getNames()}
     * @return the corresponding {@code FitSpec} subclass
     * @throws IllegalArgumentException if {@code methodName} is not registered
     */
    public static Class<? extends FitSpec> getFitSpecClass(String methodName) {
        Class<? extends FitSpec> clazz = CLASSES.get(methodName);
        if (clazz != null) {
            return clazz;
        }
        throw new IllegalArgumentException(
            String.format(
                "`methodName` should be one of: %s",
                String.join(", ", getNames())
            )
        );
    }

    // ── Fitting ─────────────────────────────────────────────────────────

    /**
     * Executes a single optimization run against the given model using the
     * supplied {@link RelaxFit} context.
     *
     * @param relaxFit configured fitter holding the experimental data
     * @param model    the model to fit
     * @return a {@link Score} summarizing the fit quality and best-fit parameters
     * @throws RuntimeException if the optimizer fails to converge
     */
    protected Score runFit(RelaxFit relaxFit, MFModelIso model) {
        double[] start = model.getStart();
        double[] lower = getLower(model);
        double[] upper = model.getUpper();

        Optional<PointValuePair> result = relaxFit.fitResidueToModel(start, lower, upper);
        if (result.isEmpty()) {
            throw new RuntimeException("Could not generate fit result.");
        }
        return relaxFit.score(result.get().getPoint(), true);
    }

    // ── TOML serialization ──────────────────────────────────────────────

    /**
     * Creates a {@link StringBuilder} pre-populated with the common
     * {@code [fit_specification]} TOML fields. Subclasses should call this
     * from {@link #toToml()} and append any additional fields.
     *
     * @return a builder containing the base TOML key-value pairs
     */
    protected StringBuilder getBaseTomlBuilder() {
        StringBuilder builder = new StringBuilder("[fit_specification]\n");
        builder.append(
            String.format(
                "method = \"%s\"%n",
                this.getClass()
                    .getSimpleName()
                    .toLowerCase()
                    .replace("fitspec", "")
            )
        );
        builder.append(String.format("tauM = %.15g%n", tauM));
        builder.append(String.format("fitTauM = %b%n", fitTauM));
        builder.append(String.format("tauMFraction = %.15g%n", tauMFraction));
        builder.append(String.format("r2Limit = %.15g%n", r2Limit));
        builder.append(String.format("fitJ = %b%n", fitJ));
        builder.append(String.format("bootstrapMode = \"%s\"%n", bootstrapMode.toString().toLowerCase()));
        builder.append(String.format("nReplicates = %d%n", nReplicates));
        builder.append(String.format("useMedian = %b%n", useMedian));
        return builder;
    }

    // ── tau_M lifecycle ─────────────────────────────────────────────────

    /**
     * Returns whether tau_M still needs to be estimated from the data.
     *
     * @return {@code true} if tau_M has not yet been set
     */
    public boolean tauMNeedsComputing() {
        return tauMNeedsComputing;
    }

    /**
     * Sets the overall tumbling correlation time and marks it as computed.
     *
     * @param tauM the correlation time in seconds; must be positive
     * @throws IllegalArgumentException if {@code tauM} is not positive
     */
    // IMPROVED: added validation
    public void setTauM(double tauM) {
        if (tauM <= 0) {
            throw new IllegalArgumentException("tauM must be positive, got: " + tauM);
        }
        this.tauM = tauM;
        this.tauMNeedsComputing = false;
    }

    /**
     * Returns the overall tumbling correlation time.
     *
     * @return tau_M in seconds
     * @throws IllegalStateException if tau_M has not been set or computed yet
     */
    public double getTauM() {
        if (tauMNeedsComputing) {
            throw new IllegalStateException("tauM has not been set!");
        }
        return tauM;
    }

    /**
     * Determines whether tau_M should be treated as a free parameter for
     * the given residue. Returns {@code true} only if global tau_M fitting
     * is enabled <em>and</em> at least one R2 value in the data exceeds
     * {@link #r2Limit}.
     *
     * @param data relaxation data for the residue under consideration
     * @return {@code true} if tau_M should be fit for this residue
     */
    public boolean fitTauM(MolDataValues<?> data) {
        return fitTauM && data.getData().stream().anyMatch(value -> value.R2 > r2Limit);
    }

    // ── Enums ───────────────────────────────────────────────────────────

    /**
     * Resampling strategy used for bootstrap error estimation.
     *
     * <ul>
     *   <li>{@link #PARAMETRIC} — resample from the fitted residual distribution</li>
     *   <li>{@link #NONPARAMETRIC} — resample observations with replacement</li>
     *   <li>{@link #BAYESIAN} — Bayesian bootstrap (Dirichlet-weighted resampling)</li>
     * </ul>
     */
    public enum BootstrapMode {
        PARAMETRIC,
        NONPARAMETRIC,
        BAYESIAN;

        /**
         * Returns a title-cased representation (e.g. "Parametric").
         *
         * @return the display name
         */
        @Override
        public String toString() {
            String name = name();
            return String.format("%c%s", name.charAt(0), name.substring(1).toLowerCase());
        }
    }

    /**
     * The chemical moiety type being analyzed, which determines the
     * appropriate spectral-density mapping, interaction constants, and
     * model-free equations.
     */
    public enum MoietyType {
        /** Backbone amide group (). */
        AMIDE("¹⁵N¹H", ""),
        /** Deuterated methyl group (). */
        DEUTERATED_METHYL("¹³C²H¹H₂", "D");

        private final String label;
        private final String modelPrefix;

        MoietyType(String label, String modelPrefix) {
            this.label = label;
            this.modelPrefix = modelPrefix;
        }

        /**
         * Returns a human-readable representation of the moiety.
         *
         * @return the display name
         */
        @Override
        public String toString() { return label; }

        public String getModelName(String baseName) { return modelPrefix + baseName; }
    }

    // ── Bootstrap ───────────────────────────────────────────────────────

    /**
     * Creates the appropriate {@link BootstrapSampler} for the configured
     * bootstrap mode.
     *
     * @param data the relaxation data to be resampled
     * @return a sampler matching {@link #bootstrapMode}
     */
    public <T extends RelaxDataValue> BootstrapSampler<T> getBootstrapSampler(MolDataValues<T> data) {
        return switch (bootstrapMode) {
            case PARAMETRIC    -> new ParametricSampler<>(data);
            case NONPARAMETRIC -> new NonparametricSampler<>(data);
            case BAYESIAN      -> new BayesianSampler<>(data);
        };
    }

    /**
     * Initializes a {@link RelaxFit} instance configured with the spectral
     * density fitting flag and the given relaxation data.
     *
     * @param key  identifier for the data set
     * @param data measured relaxation values
     * @return a ready-to-use {@link RelaxFit}
     */
    protected RelaxFit initRelaxFit(String key, MolDataValues<?> data) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setFitJ(fitJ);
        relaxFit.setRelaxData(key, data);
        return relaxFit;
    }

    // ── Statistics ───────────────────────────────────────────────────────

    /**
     * Computes central estimates and errors from bootstrap replicate parameters,
     * dispatching to the appropriate strategy based on {@link #bootstrapMode}.
     *
     * <p>For {@link BootstrapMode#PARAMETRIC}, errors are simple standard
     * deviations of the replicate distribution. For
     * {@link BootstrapMode#NONPARAMETRIC} and {@link BootstrapMode#BAYESIAN},
     * errors are derived from the covariance between bootstrap weights and
     * parameter estimates (the "infinitesimal jackknife" variance estimator).</p>
     *
     * @param parameters {@code [nParameters][nReplicates]} array of fitted values
     * @param weights    {@code [nDatapoints][nReplicates]} array of bootstrap
     *                   weights; will be {@code null} for parametric mode
     * @return a pair of (estimates, errors) arrays, each of length nParameters
     */
    protected Pair<double[], double[]> computeStatistics(double[][] parameters, double[][] weights) {
        return switch (bootstrapMode) {
            case PARAMETRIC -> computeStatisticsParametric(parameters);
            case NONPARAMETRIC, BAYESIAN -> computeStatisticsNonparametric(parameters, weights);
        };
    }

    /**
     * Returns either the mean or median of the given statistics, depending on
     * the {@link #useMedian} flag.
     *
     * @param stats descriptive statistics over a set of replicate values
     * @return the chosen central estimator
     */
    private double getAverage(DescriptiveStatistics stats) {
        return useMedian ? stats.getPercentile(50.0) : stats.getMean();
    }

    /**
     * Computes parameter estimates and standard-deviation errors from
     * parametric bootstrap replicates.
     *
     * @param parameters {@code [nParameters][nReplicates]} fitted values
     * @return pair of (estimates, errors)
     */
    private Pair<double[], double[]> computeStatisticsParametric(double[][] parameters) {
        int nParameters = parameters.length;
        double[] parameterEstimates = new double[nParameters];
        double[] parameterErrors = new double[nParameters];
        for (int k = 0; k < nParameters; k++) {
            DescriptiveStatistics stats = new DescriptiveStatistics(parameters[k]);
            parameterEstimates[k] = getAverage(stats);
            parameterErrors[k] = stats.getStandardDeviation();
        }
        return Pair.of(parameterEstimates, parameterErrors);
    }

    /**
     * Computes parameter estimates and infinitesimal-jackknife errors from
     * nonparametric or Bayesian bootstrap replicates.
     *
     * <p>The error for each parameter is computed as the square root of the
     * sum of squared covariances between the bootstrap weights and the
     * parameter values across all data points.</p>
     *
     * @param parameters {@code [nParameters][nReplicates]} fitted values
     * @param weights    {@code [nDatapoints][nReplicates]} bootstrap weights
     * @return pair of (estimates, errors)
     */
    private Pair<double[], double[]> computeStatisticsNonparametric(double[][] parameters, double[][] weights) {
        int nParameters = parameters.length;
        int nReps = parameters[0].length;
        int nDatapoints = weights.length;
        double[] parameterEstimates = new double[nParameters];
        double[] parameterErrors = new double[nParameters];

        double[] weightMeans = new double[nDatapoints];
        for (int j = 0; j < nDatapoints; j++) {
            DescriptiveStatistics weightStats = new DescriptiveStatistics(weights[j]);
            weightMeans[j] = weightStats.getMean();
        }

        for (int k = 0; k < nParameters; k++) {
            DescriptiveStatistics parameterStats = new DescriptiveStatistics(parameters[k]);
            parameterEstimates[k] = getAverage(parameterStats);
            double parameterMean = parameterStats.getMean();
            double variance = 0.0;
            for (int j = 0; j < nDatapoints; j++) {
                double weightMean = weightMeans[j];
                double covjk = 0.0;
                for (int i = 0; i < nReps; i++) {
                    covjk += (weights[j][i] - weightMean) * (parameters[k][i] - parameterMean);
                }
                covjk /= nReps;
                variance += covjk * covjk;
            }
            parameterErrors[k] = Math.sqrt(variance);
        }
        return Pair.of(parameterEstimates, parameterErrors);
    }

    // ── Order parameter construction ────────────────────────────────────

    /**
     * Constructs an {@link OrderPar} from the fit results and attaches it
     * (along with the corresponding {@link SpectralDensity}) to the atom.
     *
     * <p><strong>Side effects:</strong> This method calls
     * {@link Atom#addOrderPar(OrderParSet, OrderPar)} and
     * {@link Atom#addSpectralDensity(String, SpectralDensity)} on the atom
     * contained within {@code data}.</p>
     *
     * @param orderParSet the set to which the new order parameter belongs
     * @param data        relaxation data for the residue
     * @param key         data-set identifier
     * @param score       fit quality metrics
     * @param model       the model that was selected
     * @param parameters  best-fit parameter values
     * @param errors      parameter error estimates
     * @return the newly created {@link OrderPar}
     */
    protected OrderPar makeOrderPar(
        OrderParSet orderParSet,
        MolDataValues<?> data,
        String key,
        Score score,
        MFModelIso model,
        double[] parameters,
        double[] errors
    ) {
        ResonanceSource resSource = new ResonanceSource(data.getAtom());
        Atom atom = data.getAtom();
        List<String> names = model.getParNames();
        OrderPar orderPar = new OrderPar(
            orderParSet, resSource, score.rss, score.nValues, score.nPars, model.getName()
        );

        int nParameters = parameters.length;
        for (int k = 0; k < nParameters; k++) {
            orderPar = orderPar.set(names.get(k), parameters[k], errors[k]);
        }

        // If tauM is fixed, set it to the fixed value with zero error
        if (!model.fitTau()) {
            orderPar = orderPar.set("Tau_e", model.getTau(), 0.0);
        }

        orderPar = orderPar.set("model", (double) model.getNumber(), null);
        atom.addOrderPar(orderParSet, orderPar);

        SpectralDensity spectralDensity = new SpectralDensity(key, data.getJValues());
        atom.addSpectralDensity(key, spectralDensity);

        return orderPar;
    }

    // ── Fingerprinting ──────────────────────────────────────────────────

    /**
     * Builds a canonical, deterministic string representation of this object's
     * full state, suitable for hashing. Subclasses contribute additional fields
     * via {@link #appendSubclassState(StringBuilder)}.
     *
     * <p>The string is pipe-delimited and encodes every configuration field
     * in a fixed order. Floating-point values are encoded via
     * {@link Double#doubleToLongBits(double)} to ensure exact bit-level
     * reproducibility.</p>
     *
     * @return pipe-delimited string encoding every configuration field
     */
    protected String canonicalStateString() {
        StringBuilder sb = new StringBuilder();
        // Include concrete class so different subclasses produce different strings
        sb.append(this.getClass().getName()).append('|');

        // Core FitSpec fields in a fixed order
        sb.append("tauM=").append(Double.doubleToLongBits(tauM)).append('|');
        sb.append("tauMNeedsComputing=").append(tauMNeedsComputing).append('|');
        sb.append("fitTauM=").append(fitTauM).append('|');
        sb.append("fitJ=").append(fitJ).append('|');
        sb.append("bootstrapMode=").append(bootstrapMode == null ? "null" : bootstrapMode.name()).append('|');
        sb.append("fitExchange=").append(fitExchange).append('|');
        sb.append("tauMFraction=").append(Double.doubleToLongBits(tauMFraction)).append('|');
        sb.append("r2Limit=").append(Double.doubleToLongBits(r2Limit)).append('|');
        sb.append("nReplicates=").append(nReplicates).append('|');

        // Hook for subclasses
        appendSubclassState(sb);

        return sb.toString();
    }

    /**
     * Hook for subclasses to append their additional fields to the canonical
     * state string in a stable, deterministic order. The default implementation
     * does nothing.
     *
     * @param sb the builder to append to
     */
    protected void appendSubclassState(StringBuilder sb) {
        // default: nothing
    }

    /**
     * Computes a SHA-256 hex fingerprint of the canonical state string and
     * returns the first {@code hexLength} hex characters.
     *
     * @param hexLength number of hex characters to return (1–64 inclusive)
     * @return truncated hex digest
     * @throws IllegalArgumentException if {@code hexLength} is out of range
     */
    protected String stateFingerprintHex(int hexLength) {
        if (hexLength <= 0 || hexLength > 64) {
            throw new IllegalArgumentException("hexLength must be between 1 and 64, got: " + hexLength);
        }
        byte[] digest;
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            digest = md.digest(canonicalStateString().getBytes(StandardCharsets.UTF_8));
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("SHA-256 algorithm not available", e);
        }

        // IMPROVED: use HexFormat (Java 17+) for cleaner hex conversion
        String hex = HexFormat.of().formatHex(digest);
        return hex.substring(0, Math.min(hexLength, hex.length()));
    }

    /**
     * Produces a filename incorporating a hex fingerprint of this specification's
     * state, useful for caching or reproducibility tracking.
     *
     * <p>Example: {@code filenameWithFingerprint("result", 16, "toml")} might
     * return {@code "result-a3f8c0e912b4d567.toml"}.</p>
     *
     * @param baseName  the filename stem; must not be {@code null} or empty
     * @param hexLength number of hex characters for the fingerprint (recommended ≥ 16)
     * @param extension file extension without the leading dot; may be {@code null} or empty
     * @return the composed filename
     * @throws IllegalArgumentException if {@code baseName} is null or empty
     */
    public String filenameWithFingerprint(String baseName, int hexLength, String extension) {
        if (baseName == null || baseName.isEmpty()) {
            throw new IllegalArgumentException("baseName must not be null or empty");
        }
        String fp = stateFingerprintHex(hexLength);
        String ext = (extension == null || extension.isEmpty()) ? "" : "." + extension;
        return String.format("%s-%s%s", baseName, fp, ext);
    }

    // ── Builder ─────────────────────────────────────────────────────────

    /**
     * Abstract fluent builder for {@code FitSpec} subclasses, using the
     * curiously-recurring template pattern (CRTP) to preserve return types
     * across the inheritance hierarchy.
     *
     * <p>All fields have sensible defaults; only those that differ from the
     * defaults need to be set explicitly. Call {@link #build()} on a concrete
     * subclass builder to obtain the final {@code FitSpec} instance.</p>
     *
     * <h2>Example usage (assuming a concrete subclass {@code ConventionalFitSpec}):</h2>
     * <pre>{@code
     * FitSpec spec = new ConventionalFitSpec.Builder()
     *     .moietyType(MoietyType.AMIDE)
     *     .tauM(17.5)
     *     .fitTauM(true)
     *     .bootstrapMode(BootstrapMode.NONPARAMETRIC)
     *     .nReplicates(100)
     *     .build();
     * }</pre>
     *
     * @param <T> the concrete builder subclass type, enabling fluent method
     *            chaining without losing the concrete return type
     */
    public static abstract class Builder<T extends Builder<T>> {

        /** Default moiety: standard amide N-H backbone. */
        private static final MoietyType DEFAULT_MOIETY_TYPE = MoietyType.AMIDE;

        /** By default, tau_M is included as a free parameter. */
        private static final boolean DEFAULT_FIT_TAUM = true;

        /** Default value for fitJ. Currently always {@code true} (not yet configurable). */
        private static final boolean DEFAULT_FIT_J = true;

        /**
         * Default tauM fraction. {@link #tauMFraction}.
         */
        private static final double DEFAULT_TAUM_FRACTION = 0.25;

        /**
         * Default R2 threshold for tau_M fitting. A value of 0.0 means the
         * R2 check is effectively disabled (any positive R2 qualifies).
         */
        private static final double DEFAULT_R2_LIMIT = 0.0;

        /** Default bootstrap resampling strategy. */
        private static final BootstrapMode DEFAULT_BOOTSTRAP_MODE = BootstrapMode.PARAMETRIC;

        /**
         * Default value for fitExchange. Currently always {@code false}
         * (not yet configurable).
         */
        private static final boolean DEFAULT_FIT_EXCHANGE = false;

        /**
         * Default number of bootstrap replicates.
         *
         * <p>25 is chosen because it is smaller than 27, the maximum number of
         * unique resamples possible for a 2-field amide dataset under
         * nonparametric bootstrapping. This avoids an unexpected error when
         * the user runs with the default set. In general, a great number of
         * replicates should be generated, if possible.</p>
         */
        private static final int DEFAULT_N_REPLICATES = 25;

        /** By default, use the mean (not median) as the central estimator. */
        private static final boolean DEFAULT_USE_MEDIAN = false;

        // ── Mutable state ───────────────────────────────────────────────

        protected MoietyType moietyType = DEFAULT_MOIETY_TYPE;
        protected Optional<Double> tauM = Optional.empty();
        protected boolean fitTauM = DEFAULT_FIT_TAUM;
        protected boolean fitJ = DEFAULT_FIT_J;
        protected BootstrapMode bootstrapMode = DEFAULT_BOOTSTRAP_MODE;
        protected boolean fitExchange = DEFAULT_FIT_EXCHANGE;
        protected double tauMFraction = DEFAULT_TAUM_FRACTION;
        protected double r2Limit = DEFAULT_R2_LIMIT;
        protected int nReplicates = DEFAULT_N_REPLICATES;
        protected boolean useMedian = DEFAULT_USE_MEDIAN;

        // ── Default-value accessors ─────────────────────────────────────

        /** Returns the default moiety type ({@link MoietyType#AMIDE}). */
        public static MoietyType getDefaultMoietyType() { return DEFAULT_MOIETY_TYPE; }

        /** Returns the default value for the fit-tau_M flag ({@code true}). */
        public static boolean getDefaultFitTauM() { return DEFAULT_FIT_TAUM; }

        /** Returns the default tau_M fraction tolerance (0.25). */
        public static double getDefaultTauMFraction() { return DEFAULT_TAUM_FRACTION; }

        /** Returns the default R2 limit (0.0). */
        public static double getDefaultR2Limit() { return DEFAULT_R2_LIMIT; }

        /** Returns the default number of bootstrap replicates (25). */
        public static int getDefaultNReplicates() { return DEFAULT_N_REPLICATES; }

        /** Returns the default bootstrap mode ({@link BootstrapMode#PARAMETRIC}). */
        public static BootstrapMode getDefaultBootstrapMode() { return DEFAULT_BOOTSTRAP_MODE; }

        /** Returns the default useMedian flag ({@code false}). */
        public static boolean getDefaultUseMedian() { return DEFAULT_USE_MEDIAN; }

        // ── CRTP self-type helper ───────────────────────────────────────

        /**
         * Returns {@code this} cast to the concrete builder type {@code T}.
         * This unchecked cast is safe as long as subclasses correctly
         * parameterize {@code T} with their own type.
         *
         * @return this builder as type {@code T}
         */
        @SuppressWarnings("unchecked")
        protected T self() {
            return (T) this;
        }

        // ── Fluent setters ──────────────────────────────────────────────

        /**
         * Sets the chemical moiety type.
         *
         * @param moietyType the moiety type (e.g. {@link MoietyType#AMIDE});
         *                   must not be {@code null}
         * @return this builder
         * @throws NullPointerException if {@code moietyType} is {@code null}
         */
        // IMPROVED: added null check
        public T moietyType(MoietyType moietyType) {
            this.moietyType = Objects.requireNonNull(moietyType, "moietyType must not be null");
            return self();
        }

        /**
         * Sets the overall tumbling correlation time to a known value.
         *
         * @param tauM correlation time in seconds; must be positive
         * @return this builder
         * @throws IllegalArgumentException if {@code tauM} is not positive
         */
        public T tauM(double tauM) {
            if (tauM <= 0) {
                throw new IllegalArgumentException("tauM must be positive, got: " + tauM);
            }
            this.tauM = Optional.of(tauM);
            return self();
        }

        /**
         * Clears the tau_M value, indicating that it should be estimated
         * from the data.
         *
         * @return this builder
         */
        public T tauMUnset() {
            this.tauM = Optional.empty();
            return self();
        }

        /**
         * Sets whether tau_M should be included as a free fit parameter.
         *
         * @param fitTauM {@code true} to fit tau_M per residue
         * @return this builder
         */
        public T fitTauM(boolean fitTauM) {
            this.fitTauM = fitTauM;
            return self();
        }

        /**
         * Sets whether spectral-density J values should be fit directly.
         *
         * <p><strong>Note:</strong> This feature is not yet implemented.
         * The supplied value is currently ignored and {@code fitJ} is
         * always set to {@code true}.</p>
         *
         * @param fitJ desired setting (currently ignored)
         * @return this builder
         */
        public T fitJ(boolean fitJ) {
            // TODO: implement — currently hard-coded to true
            this.fitJ = DEFAULT_FIT_J;
            return self();
        }

        /**
         * Sets the bootstrap resampling mode used for error estimation.
         *
         * @param bootstrapMode the resampling strategy; must not be {@code null}
         * @return this builder
         * @throws NullPointerException if {@code bootstrapMode} is {@code null}
         */
        // IMPROVED: added null check
        public T bootstrapMode(BootstrapMode bootstrapMode) {
            this.bootstrapMode = Objects.requireNonNull(bootstrapMode, "bootstrapMode must not be null");
            return self();
        }

        /**
         * Sets whether chemical-exchange contributions should be fit.
         *
         * <p><strong>Note:</strong> This feature is not yet implemented.
         * The supplied value is currently ignored and {@code fitExchange}
         * is always set to {@code false}.</p>
         *
         * @param fitExchange desired setting (currently ignored)
         * @return this builder
         */
        public T fitExchange(boolean fitExchange) {
            // TODO: implement — currently hard-coded to false
            this.fitExchange = DEFAULT_FIT_EXCHANGE;
            return self();
        }

        /**
         * Sets the tauM fractional. See {@link #tauMFraction}.
         *
         * @param tauMFraction value in the range [0, 1]
         * @return this builder
         * @throws IllegalArgumentException if {@code tauMFraction} is not in the required range.
         */
        public T tauMFraction(double tauMFraction) {
            if (tauMFraction < 0 || tauMFraction > 1) {
                throw new IllegalArgumentException("tauMFraction must be between 0 and 1, got: " + tauMFraction);
            }
            this.tauMFraction = tauMFraction;
            return self();
        }

        /**
         * Sets the R2 threshold above which tau_M fitting is attempted.
         *
         * @param r2Limit the R2 cutoff value (s⁻¹); must be non-negative
         * @return this builder
         * @throws IllegalArgumentException if {@code r2Limit} is negative
         */
        public T r2Limit(double r2Limit) {
            if (r2Limit < 0) {
                throw new IllegalArgumentException("r2Limit must be non-negative, got: " + r2Limit);
            }
            this.r2Limit = r2Limit;
            return self();
        }

        /**
         * Sets the number of bootstrap replicates for error estimation.
         *
         * @param nReplicates number of replicates; must be positive
         * @return this builder
         * @throws IllegalArgumentException if {@code nReplicates} is not positive
         */
        public T nReplicates(int nReplicates) {
            if (nReplicates <= 0) {
                throw new IllegalArgumentException("nReplicates must be positive, got: " + nReplicates);
            }
            this.nReplicates = nReplicates;
            return self();
        }

        /**
         * Sets whether the median (rather than the mean) should be used as
         * the central estimator when aggregating bootstrap replicates.
         *
         * <p><strong>Note:</strong> This setting is only meaningful for
         * {@link BootstrapMode#NONPARAMETRIC} and {@link BootstrapMode#BAYESIAN}
         * modes. For {@link BootstrapMode#PARAMETRIC}, it is ignored by the
         * Conventional fitting strategy.</p>
         *
         * @param useMedian {@code true} to use the median
         * @return this builder
         */
        public T useMedian(boolean useMedian) {
            this.useMedian = useMedian;
            return self();
        }

        /**
         * Validates cross-field constraints. Called by {@link #build()} before
         * constructing the {@code FitSpec} instance.
         *
         * <p>Subclass builders may override to add additional checks but should
         * call {@code super.validate()} first.</p>
         *
         * @throws IllegalStateException if the configuration is invalid
         */
        protected void validate() {
            Objects.requireNonNull(moietyType, "moietyType must not be null");
            Objects.requireNonNull(bootstrapMode, "bootstrapMode must not be null");
            Objects.requireNonNull(tauM, "tauM Optional must not be null");
            tauM.ifPresent(v -> {
                if (v <= 0) {
                    throw new IllegalStateException("tauM must be positive when set, got: " + v);
                }
            });
            if (nReplicates <= 0) {
                throw new IllegalStateException("nReplicates must be positive, got: " + nReplicates);
            }
            if (tauMFraction < 0 || tauMFraction > 1) {
                throw new IllegalStateException("tauMFraction must be between 0 and 1, got: " + tauMFraction);
            }
            if (r2Limit < 0) {
                throw new IllegalStateException("r2Limit must be non-negative, got: " + r2Limit);
            }
        }

        /**
         * Validates the configuration and constructs the concrete
         * {@code FitSpec} instance.
         *
         * <p>Implementations should call {@link #validate()} before
         * constructing the instance:</p>
         * <pre>{@code
         * @Override
         * public MyFitSpec build() {
         *     validate();
         *     return new MyFitSpec(this);
         * }
         * }</pre>
         *
         * @return the fully configured {@code FitSpec}
         * @throws IllegalStateException if the configuration is invalid
         */
        public abstract FitSpec build();
    }
}
