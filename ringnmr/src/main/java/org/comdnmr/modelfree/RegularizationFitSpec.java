package org.comdnmr.modelfree;

import java.util.Arrays;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;

/**
 * Regularized model-free fitting strategy using the extended
 * ({@code 2sf}) model exclusively.
 *
 * <p>Unlike {@link BaggingFitSpec}, which selects the best model per replicate
 * from a candidate set, this strategy always fits the {@code 2sf} model
 * (parameters: {@code [τm?,] S²f, τf, S²s, τs}) with
 * regularization terms that penalize deviations of the order parameters and
 * correlation times from physically motivated null values. The regularization
 * discourages overfitting without the need for AICc-based model selection.</p>
 *
 * <p>The three regularization strengths — {@link #lambdaS}, {@link #lambdaTauF},
 * and {@link #lambdaTauS} — are passed through to {@link RelaxFit} via
 * {@link RelaxFit#setLambdaS(double)}, {@link RelaxFit#setLambdaTauF(double)},
 * and {@link RelaxFit#setLambdaTauS(double)} respectively. All lower bounds for
 * the model parameters are set to zero (see {@link #getLower(MFModelIso)}).</p>
 *
 * <p>Bootstrap resampling and statistics follow the same procedure as
 * {@link BaggingFitSpec}: for each replicate, the data are resampled via the
 * configured {@link FitSpec.BootstrapMode}, the 2sf model is fit, and the
 * resulting parameters are post-processed via
 * {@link #processParamsAfterFit(double[], boolean)} before being aggregated
 * across replicates by
 * {@link FitSpec#computeStatistics(double[][], double[][])}.</p>
 *
 * <h2>Example usage:</h2>
 * <pre>{@code
 * FitSpec spec = new RegularizationFitSpec.Builder()
 *     .tauM(17.5)
 *     .bootstrapMode(BootstrapMode.NONPARAMETRIC)
 *     .nReplicates(200)
 *     .lambdaS(0.5)
 *     .lambdaTauF(0.1)
 *     .lambdaTauS(0.2)
 *     .build();
 *
 * ModelFitResult result = spec.fit(key, data, orderParSetMap);
 * }</pre>
 *
 * @see BaggingFitSpec
 * @see ConventionalFitSpec
 * @see FitSpec.BootstrapMode
 */
public class RegularizationFitSpec extends FitSpec {

    /** Key used when registering results in the {@link OrderParSet} map. */
    private static final String KEY = "REGULARIZATION";

    /**
     * Order-parameter threshold above which the corresponding motion is
     * considered absent (S² ≈ 1). When both S² values exceed this limit,
     * the result is mapped to "model 0" (no local motions).
     */
    static final double S2_THOLD = 0.99;

    /**
     * Correlation-time threshold (ns) below which the corresponding motion is
     * treated as instantaneous (τ → 0). Motions with τ less than this limit
     * are assigned τ = 0.
     */
    static final double TAU_THOLD = 1.0e-3;  // 1 ps

    /**
     * Regularization strength for the order parameters (S²f, S²s).
     * Higher values push S² estimates towards 1 (rigid limit).
     */
    private final double lambdaS;

    /**
     * Regularization strength for the fast correlation time τ_f.
     * Higher values push τf estimates towards 0.
     */
    private final double lambdaTauF;

    /**
     * Regularization strength for the slow correlation time τs.
     * Higher values suppress spurious slow-motion contributions.
     */
    private final double lambdaTauS;

    /**
     * Builder for {@link RegularizationFitSpec}.
     *
     * <p>All base configuration (tau_M, bootstrap mode, number of replicates,
     * etc.) is inherited from {@link FitSpec.Builder}. The three
     * regularization strengths have sensible defaults and may be set
     * individually.</p>
     *
     * <h2>Example:</h2>
     * <pre>{@code
     * FitSpec spec = new RegularizationFitSpec.Builder()
     *     .tauM(17.5)
     *     .bootstrapMode(BootstrapMode.NONPARAMETRIC)
     *     .nReplicates(200)
     *     .lambdaS(0.5)
     *     .lambdaTauF(0.1)
     *     .lambdaTauS(0.2)
     *     .build();
     * }</pre>
     */
    public static class Builder extends FitSpec.Builder<Builder> {

        private static final double DEFAULT_LAMBDA_S    = 0.5;
        private static final double DEFAULT_LAMBDA_TAUF = 0.1;
        private static final double DEFAULT_LAMBDA_TAUS = 0.2;

        private double lambdaS    = DEFAULT_LAMBDA_S;
        private double lambdaTauF = DEFAULT_LAMBDA_TAUF;
        private double lambdaTauS = DEFAULT_LAMBDA_TAUS;

        /** Returns the default regularization strength for S² (0.5). */
        public static double getDefaultLambdaS()    { return DEFAULT_LAMBDA_S; }

        /** Returns the default regularization strength for τ_f (0.1). */
        public static double getDefaultLambdaTauF() { return DEFAULT_LAMBDA_TAUF; }

        /** Returns the default regularization strength for τ_s (0.2). */
        public static double getDefaultLambdaTauS() { return DEFAULT_LAMBDA_TAUS; }

        private void validateLambda(String name, double value) {
            if (value < 0.0) {
                throw new IllegalArgumentException(
                    String.format("%s must be >= 0.0, got: %s", name, value)
                );
            }
        }

        /**
         * Sets the regularization strength for the order parameters.
         *
         * @param lambdaS regularization weight; must be &ge; 0
         * @return this builder
         * @throws IllegalArgumentException if {@code lambdaS} is negative
         */
        public Builder lambdaS(double lambdaS) {
            validateLambda("lambdaS", lambdaS);
            this.lambdaS = lambdaS;
            return this;
        }

        /**
         * Sets the regularization strength for the fast correlation time τ_f.
         *
         * @param lambdaTauF regularization weight; must be &ge; 0
         * @return this builder
         * @throws IllegalArgumentException if {@code lambdaTauF} is negative
         */
        public Builder lambdaTauF(double lambdaTauF) {
            validateLambda("lambdaTauF", lambdaTauF);
            this.lambdaTauF = lambdaTauF;
            return this;
        }

        /**
         * Sets the regularization strength for the slow correlation time τ_s.
         *
         * @param lambdaTauS regularization weight; must be &ge; 0
         * @return this builder
         * @throws IllegalArgumentException if {@code lambdaTauS} is negative
         */
        public Builder lambdaTauS(double lambdaTauS) {
            validateLambda("lambdaTauS", lambdaTauS);
            this.lambdaTauS = lambdaTauS;
            return this;
        }

        /**
         * Validates the configuration and constructs a new
         * {@link RegularizationFitSpec}.
         *
         * @return the configured {@code RegularizationFitSpec}
         * @throws IllegalStateException if the configuration is invalid
         */
        @Override
        public RegularizationFitSpec build() {
            validate();
            return new RegularizationFitSpec(this);
        }
    }

    /**
     * Constructs a new {@code RegularizationFitSpec} from the given builder.
     *
     * @param builder the builder containing all configuration values
     */
    protected RegularizationFitSpec(Builder builder) {
        super(builder);
        this.lambdaS    = builder.lambdaS;
        this.lambdaTauF = builder.lambdaTauF;
        this.lambdaTauS = builder.lambdaTauS;
    }

    /** Returns the order-parameter regularization strength. */
    double getLambdaS()    { return lambdaS; }

    /** Returns the fast-correlation-time regularization strength. */
    double getLambdaTauF() { return lambdaTauF; }

    /** Returns the slow-correlation-time regularization strength. */
    double getLambdaTauS() { return lambdaTauS; }

    /**
     * {@inheritDoc}
     *
     * <p>All lower bounds are set to zero, allowing the regularization terms
     * to act as soft constraints rather than relying on hard box bounds.</p>
     */
    @Override
    protected double[] getLower(MFModelIso model) {
        return new double[model.getNPars()];
    }

    /**
     * {@inheritDoc}
     *
     * <p>Appends {@code lambdaS}, {@code lambdaTauF}, and {@code lambdaTauS}
     * to the canonical state string.</p>
     */
    @Override
    protected void appendSubclassState(StringBuilder sb) {
        sb.append("lambdaS=").append(Double.doubleToLongBits(lambdaS)).append('|');
        sb.append("lambdaTauF=").append(Double.doubleToLongBits(lambdaTauF)).append('|');
        sb.append("lambdaTauS=").append(Double.doubleToLongBits(lambdaTauS)).append('|');
    }

    @Override
    public String toToml() {
        StringBuilder builder = getBaseTomlBuilder();
        builder.append(String.format("lambdaS = %s%n", lambdaS));
        builder.append(String.format("lambdaTauF = %s%n", lambdaTauF));
        builder.append(String.format("lambdaTauS = %s", lambdaTauS));
        return builder.toString();
    }

    // ── RelaxFit initialization ───────────────────────────────────────────────

    /**
     * {@inheritDoc}
     *
     * <p>In addition to the base initialization, enables lambda regularization
     * on the {@link RelaxFit} instance and sets the three regularization
     * strengths.</p>
     */
    @Override
    protected RelaxFit initRelaxFit(String key, MolDataValues<? extends RelaxDataValue> data) {
        RelaxFit relaxFit = super.initRelaxFit(key, data);
        relaxFit.setUseLambda(true);
        relaxFit.setLambdaS(getLambdaS());
        relaxFit.setLambdaTauF(getLambdaTauF());
        relaxFit.setLambdaTauS(getLambdaTauS());
        return relaxFit;
    }

   /**
    * Processes output of 2sf fit.
    *
    * <ul>
    *   <li><em>No local motions</em> (both S² above {@code S2_THOLD}):
    *       S²f = S²s = 1, τf = τs = 0.</li>
    *   <li><em>One fast motion</em> (τ below {@code TAU_THOLD}):
    *       S²f = S², τf = 0, S²s = 1, τs = 0.</li>
    *   <li><em>One fast motion with non-zero τ_f</em>
    *       (τ in (TAU_THOLD, SLOW_LIMIT]):
    *       S²f = S², τf = τ, S²s = 1, τs = 0.</li>
    *   <li><em>One slow motion</em> (τ above {@code SLOW_LIMIT}):
    *       S²f = 1, τf = 0, S²s = S², τs = τ.</li>
    *   <li><em>Two motions, one instantaneous</em>:
    *       The instantaneous timescale is mapped to τf = 0.</li>
    *   <li><em>Two independent motions</em>:
    *       Sorted so that τf &lt; τs.</li>
    * </ul>
    *
    * @param params parameter values from the optimizer
    * @return the processed parameter array
    */
    protected double[] processParamsAfterFit(double[] params, boolean fitTau) {
       return params;

       // FIXME: was using this for amide data. It is incompatible for deuterium
       // int start = fitTau ? 1 : 0;
       // double s1   = params[start];
       // double tau1 = params[start + 1];
       // double s2   = params[start + 2];
       // double tau2 = params[start + 3];

       // double sf2, tauf, ss2, taus;

       // // No local motions — both order parameters are effectively 1
       // if (s1 > S2_THOLD && s2 > S2_THOLD) {
       //     sf2 = ss2 = 1.0;
       //     tauf = taus = 0.0;
       // }

       // // One local motion — the other degree of freedom is frozen out
       // else if (s1 > S2_THOLD || s2 > S2_THOLD) {
       //     double s   = (s1 > S2_THOLD) ? s2 : s1;
       //     double tau = (s1 > S2_THOLD) ? tau2 : tau1;
       //     if (tau < TAU_THOLD) {
       //         // Fast motion with effectively zero correlation time (Model 1)
       //         sf2 = s;   tauf = 0.0; ss2 = 1.0; taus = 0.0;
       //     } else if (tau < MFModelIso2sf.SLOW_LIMIT) {
       //         // Resolvable fast motion (Model 1f)
       //         sf2 = s;   tauf = tau; ss2 = 1.0; taus = 0.0;
       //     } else {
       //         // Slow motion (Model 1s)
       //         sf2 = 1.0; tauf = 0.0; ss2 = s;   taus = tau;
       //     }
       // }

       // // Two local motions
       // else {
       //     if (tau1 < TAU_THOLD) {
       //         // tau1 is instantaneous: assign it to the fast slot (Model 2s)
       //         sf2 = s1; tauf = 0.0; ss2 = s2; taus = tau2;
       //     } else if (tau2 < TAU_THOLD) {
       //         // tau2 is instantaneous: assign it to the fast slot (Model 2s)
       //         sf2 = s2; tauf = 0.0; ss2 = s1; taus = tau1;
       //     } else {
       //         // Both timescales are resolvable: sort so that tauf < taus (Model 2sf)
       //         if (tau1 < tau2) {
       //             sf2 = s1; tauf = tau1; ss2 = s2; taus = tau2;
       //         } else {
       //             sf2 = s2; tauf = tau2; ss2 = s1; taus = tau1;
       //         }
       //     }
       // }

       // params[start]     = sf2;
       // params[start + 1] = tauf;
       // params[start + 2] = ss2;
       // params[start + 3] = taus;

       // return params;
    }

    /**
     * Performs a regularized model-free fit for a single residue.
     *
     * <p>The algorithm proceeds as follows:</p>
     * <ol>
     *   <li>For each of the {@code nReplicates} bootstrap replicates:
     *     <ol type="a">
     *       <li>Draw a bootstrap sample from {@code data}.</li>
     *       <li>Fit the {@code 2sf} model to the sample with regularization.</li>
     *       <li>Post-process the parameters via
     *           {@link #processParamsAfterFit(double[], boolean)}.</li>
     *       <li>Record the parameters and bootstrap weights.</li>
     *     </ol>
     *   </li>
     *   <li>Aggregate parameters and weights across replicates using
     *       {@link #computeStatistics(double[][], double[][])} to obtain
     *       final estimates and errors.</li>
     *   <li>Construct and register an {@link OrderPar} under the key
     *       {@code "REGULARIZATION"}.</li>
     * </ol>
     *
     * <p><strong>Note on Score:</strong> The {@link Score} stored in the
     * returned {@link ModelFitResult} and passed to {@link #makeOrderPar} is
     * taken from the first replicate. A meaningful aggregate score has not yet
     * been defined.</p>
     *
     * @param key            identifier for the relaxation data set
     * @param data           measured relaxation values for the residue
     * @param orderParSetMap mutable map to which the resulting
     *                       {@link OrderParSet} entry will be added
     * @return the {@link ModelFitResult} with aggregated parameters and errors
     * @throws IllegalStateException if tau_M has not been set
     */
    @Override
    public ModelFitResult fit(String key, MolDataValues<?> data, Map<String, OrderParSet> orderParSetMap) {
        RelaxFit relaxFit = initRelaxFit(key, data);
        MFModelIso2sf model = (MFModelIso2sf) getModel("2sf", data);
        data.setTestModel(model);

        int nParameters = model.getNPars();
        int nWeights = data.getNSpectralDensities();
        double[][] parameters = new double[nParameters][nReplicates];
        double[][] weights = new double[nWeights][nReplicates];
        BootstrapSampler<? extends RelaxDataValue> sampler = getBootstrapSampler(data);

        Score[] scores = new Score[nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues<? extends RelaxDataValue> replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            scores[i] = runFit(relaxFit, model);
            double[] replicateParameters = processParamsAfterFit(scores[i].getPars(), model.fitTau());
            double[] replicateWeights = replicateData.getWeights();
            for (int k = 0; k < nParameters; k++) parameters[k][i] = replicateParameters[k];
            for (int j = 0; j < nWeights; j++) weights[j][i] = replicateWeights[j];
        }

        Pair<double[], double[]> parameterEstimates = computeStatistics(parameters, weights);
        double[] fitParameters = parameterEstimates.getLeft();
        System.out.printf("fitParameters:%n%s%n", Arrays.toString(fitParameters));
        double[] fitErrors = parameterEstimates.getRight();
        System.out.printf("fitErrors:%n%s%n", Arrays.toString(fitErrors));

        orderParSetMap.computeIfAbsent(KEY, ky -> new OrderParSet(ky));
        // FIXME: the Score used here (scores[0]) is from the first replicate.
        // For bootstrap fitting, there is no single score over the original
        // data; this should be revisited.
        OrderPar orderPar = makeOrderPar(
            orderParSetMap.get(KEY),
            sampler.getOriginalData(),
            key,
            scores[0],
            model,
            fitParameters,
            fitErrors
        );

        return new ModelFitResult(orderPar, parameters, null);
    }
}
