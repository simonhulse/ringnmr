package org.comdnmr.modelfree;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.Pair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;

/**
 * Bootstrap-aggregation (bagging) model-free fitting strategy.
 *
 * <p>For each bootstrap replicate, every candidate {@link MFModelIso} model is
 * fitted to the resampled data and the model with the lowest AICc is selected.
 * The final parameter estimates and errors are derived by aggregating the
 * optimal model parameters across all replicates (see
 * {@link FitSpec#computeStatistics(double[][], double[][])}).</p>
 *
 * <p>Because model selection happens independently per replicate, the aggregate
 * result is implicitly model-averaged: replicates may vote for different models,
 * and the per-parameter distributions naturally reflect that uncertainty. After
 * selection, each replicate's raw parameters are converted to the extended
 * model-free ({@code 2sf}) representation via
 * {@link MFModelIso#getStandardPars(double[])} before being stored.</p>
 *
 * <p>The set of candidate models is determined by the moiety type: AMIDE uses
 * models 1, 1f, 1s, 2s, 2sf; DEUTERATED_METHYL additionally includes 1sf.
 * Other common settings (tau_M, bootstrap mode, number of replicates, etc.)
 * are inherited from {@link FitSpec}.</p>
 *
 * <h2>Example usage:</h2>
 * <pre>{@code
 * FitSpec spec = new BaggingFitSpec.Builder()
 *     .tauM(17.5)
 *     .bootstrapMode(BootstrapMode.NONPARAMETRIC)
 *     .nReplicates(200)
 *     .build();
 *
 * ModelFitResult result = spec.fit(key, data, orderParSetMap);
 * }</pre>
 *
 * @see ConventionalFitSpec
 * @see RegularizationFitSpec
 */
public class BaggingFitSpec extends FitSpec {

    /** Key used when registering results in the {@link OrderParSet} map. */
    private static final String KEY = "BAGGING";

    // ── Builder ─────────────────────────────────────────────────────────────

    /**
     * Fluent builder for {@link BaggingFitSpec}.
     *
     * <p>All configuration is inherited from {@link FitSpec.Builder} (tau_M,
     * bootstrap mode, number of replicates, etc.). The model list is derived
     * automatically from the moiety type. No bagging-specific parameters are
     * required.</p>
     */
    public static class Builder extends FitSpec.Builder<Builder> {

        /**
         * Validates the configuration and constructs a new {@link BaggingFitSpec}.
         *
         * @return the configured {@code BaggingFitSpec}
         * @throws IllegalStateException if the configuration is invalid
         */
        @Override
        public BaggingFitSpec build() {
            validate();
            return new BaggingFitSpec(this);
        }
    }

    // ── Constructor ──────────────────────────────────────────────────────────

    /**
     * Constructs a new {@code BaggingFitSpec} from the given builder.
     *
     * @param builder the builder containing all configuration values
     */
    protected BaggingFitSpec(Builder builder) { super(builder); }

    // ── Model selection ───────────────────────────────────────────────────────

    public List<String> getModelNames() {
        List<String> names = new ArrayList<>(Arrays.asList("1", "1f", "1s", "2s", "2sf"));
        if (moietyType == MoietyType.DEUTERATED_METHYL) {
            names.add("1sf");
        }
        return names;
    }

    // FIXME: repeated with ConventionalFitSpec
    public List<MFModelIso> getModels(MolDataValues<? extends RelaxDataValue> data) {
        return getModelNames()
            .stream()
            .map(name -> getModel(name, data))
            .collect(Collectors.toList());
    }

    @Override
    public String toToml() { return getBaseTomlBuilder().toString(); }

    // ── fit ──────────────────────────────────────────────────────────────────

    /**
     * Performs a bagging model-free fit for a single residue.
     *
     * <p>The algorithm proceeds as follows:</p>
     * <ol>
     *   <li>For each of the {@code nReplicates} bootstrap replicates:
     *     <ol type="a">
     *       <li>Draw a bootstrap sample from {@code data}.</li>
     *       <li>Fit every candidate model to the sample.</li>
     *       <li>Select the model with the lowest AICc.</li>
     *       <li>Convert the optimal parameters to the {@code 2sf}
     *           representation using
     *           {@link MFModelIso#getStandardPars(double[])}.</li>
     *       <li>Record the parameters and bootstrap weights.</li>
     *     </ol>
     *   </li>
     *   <li>Aggregate parameters and weights across replicates using
     *       {@link #computeStatistics(double[][], double[][])} to obtain
     *       final estimates and errors.</li>
     *   <li>Construct and register an {@link OrderPar} under the key
     *       {@code "BAGGING"}.</li>
     * </ol>
     *
     * <p><strong>Note on Score:</strong> The {@link Score} stored in the
     * returned {@link ModelFitResult} and passed to {@link #makeOrderPar} is
     * taken from the first replicate's best fit. This is a known limitation —
     * a meaningful aggregate score for bagging has not yet been defined.</p>
     *
     * @param key            identifier for the relaxation data set
     * @param data           measured relaxation values for the residue
     * @param orderParSetMap mutable map to which the resulting
     *                       {@link OrderParSet} entry will be added
     * @return the {@link ModelFitResult} with aggregated parameters and errors
     * @throws AssertionError    if no model could be fit for a replicate
     * @throws IllegalStateException if tau_M has not been set
     */
    @Override
    public ModelFitResult fit(String key, MolDataValues<?> data, Map<String, OrderParSet> orderParSetMap) {
        RelaxFit relaxFit = initRelaxFit(key, data);

        List<MFModelIso> models = getModels(data);
        MFModelIso2sf model2sf = (MFModelIso2sf) getModel("2sf", data);

        // nParameters is 5 when tauM is free, 4 when it is fixed
        int nParameters = model2sf.getNPars();
        int nWeights = data.getNValues();
        double[][] parameters = new double[nParameters][nReplicates];
        double[][] weights = new double[nWeights][nReplicates];
        BootstrapSampler<? extends RelaxDataValue> sampler = getBootstrapSampler(data);

        Score[] bestScores = new Score[nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues<? extends RelaxDataValue> replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);

            Optional<Pair<Score, MFModelIso>> bestScoreModel = Optional.empty();
            for (MFModelIso model : models) {
                data.setTestModel(model);
                Score score = runFit(relaxFit, model);
                if (bestScoreModel.isEmpty() || score.aicc().get() < bestScoreModel.get().getLeft().aicc().get()) {
                    bestScoreModel = Optional.of(Pair.of(score, model));
                }
            }

            if (bestScoreModel.isEmpty()) {
                throw new AssertionError("`bestScore` should not be empty here!");
            }

            bestScores[i] = bestScoreModel.get().getLeft();
            MFModelIso bestModel = bestScoreModel.get().getRight();
            // Convert the winning model's parameters to the canonical 2sf
            // representation so that all replicates are on a common scale.
            double[] replicateParameters = bestModel.getStandardPars(bestScores[i].getPars());
            double[] replicateWeights = replicateData.getWeights();
            for (int k = 0; k < nParameters; k++) parameters[k][i] = replicateParameters[k];
            for (int j = 0; j < nWeights; j++) weights[j][i] = replicateWeights[j];
        }

        Pair<double[], double[]> parameterEstimates = computeStatistics(parameters, weights);
        double[] fitParameters = parameterEstimates.getLeft();
        double[] fitErrors = parameterEstimates.getRight();

        orderParSetMap.computeIfAbsent(KEY, ky -> new OrderParSet(ky));
        // FIXME: the Score used here (bestScores[0]) is from the first
        // replicate's best fit. For bagging, there is no single score
        // over the original data; this should be revisited.
        OrderPar orderPar = makeOrderPar(
            orderParSetMap.get(KEY),
            sampler.getOriginalData(),
            key,
            bestScores[0],
            model2sf,
            fitParameters,
            fitErrors
        );

        return new ModelFitResult(orderPar, parameters, null);
    }
}
