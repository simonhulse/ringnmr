package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import org.comdnmr.modelfree.models.MFModelIso;

import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;

/**
 * Conventional (non-aggregated) model-free fitting strategy.
 *
 * <p>This strategy fits each candidate {@link MFModelIso} independently to the
 * original data and selects the best model using the small sample-corrected Akaike
 * Information Criterion (AICc). Parameter errors are estimated via
 * bootstrap resampling (parametric, nonparametric, or Bayesian depending
 * on the configured {@link BootstrapMode}), but the final parameter
 * <em>values</em> come from the fit to the original data.</p>
 *
 * <p>By default, all registered {@link MFModelIso} models are tested. The
 * list can be customized via {@link Builder#modelNames(List)}. For
 * {@link MoietyType#AMIDE}, the {@code "1sf"} model is automatically
 * excluded.</p>
 *
 * <h2>Example usage:</h2>
 * <pre>{@code
 * FitSpec spec = new ConventionalFitSpec.Builder()
 *     .tauM(17.5)
 *     .bootstrapMode(BootstrapMode.NONPARAMETRIC)
 *     .nReplicates(100)
 *     .modelNames(List.of("1", "1f", "2s", "2sf"))
 *     .build();
 *
 * ModelFitResult result = spec.fit(key, data, orderParSetMap);
 * }</pre>
 *
 * @see BaggingFitSpec
 * @see RegularizationFitSpec
 */
public class ConventionalFitSpec extends FitSpec {

    /** Prefix used when constructing {@link OrderParSet} keys for this strategy. */
    private static final String KEY = "CONVENTIONAL";

    /**
     * Names of the {@link MFModelIso} models to evaluate during fitting.
     * Determined at construction time based on the builder configuration
     * and the {@link MoietyType}.
     */
    protected final List<String> modelNames;

    /**
     * Fluent builder for {@link ConventionalFitSpec}.
     *
     * <p>Inherits all configuration from {@link FitSpec.Builder} and adds
     * an optional {@link #modelNames(List)} setting to restrict which models
     * are evaluated.</p>
     */
    public static class Builder extends FitSpec.Builder<Builder> {

        /** The complete set of model names. */
        private static final Set<String> ALL_MODEL_NAMES =
            Set.of("1", "1f", "1s", "1sf", "2s", "2sf");

        /** Model names that are only valid for {@link MoietyType#DEUTERATED_METHYL}. */
        private static final Set<String> DEUTERATED_METHYL_ONLY =
            Set.of("1sf");

        /** Optional explicit list of model names to evaluate. */
        protected Optional<List<String>> modelNames = Optional.empty();

        /**
         * Validates the configuration and constructs a new
         * {@link ConventionalFitSpec}.
         *
         * @return the configured {@code ConventionalFitSpec}
         * @throws IllegalStateException if the configuration is invalid,
         *         including invalid model names for the configured moiety type
         */
        public ConventionalFitSpec build() {
            validate();
            return new ConventionalFitSpec(this);
        }

        /**
         * {@inheritDoc}
         *
         * <p>Additionally validates that any explicitly provided model names
         * are recognized and compatible with the configured
         * {@link MoietyType}.</p>
         */
        @Override
        protected void validate() {
            super.validate();
            modelNames.ifPresent(names -> {
                if (names.isEmpty()) {
                    throw new IllegalStateException("modelNames must not be empty");
                }
                for (String name : names) {
                    if (!ALL_MODEL_NAMES.contains(name)) {
                        throw new IllegalArgumentException(
                            String.format(
                                "Unknown model name \"%s\"; valid names are: %s",
                                name,
                                ALL_MODEL_NAMES
                            )
                        );
                    }
                    if (DEUTERATED_METHYL_ONLY.contains(name) && moietyType != MoietyType.DEUTERATED_METHYL) {
                        throw new IllegalArgumentException(
                            String.format(
                                "Model \"%s\" is only valid for %s, but moiety type is %s",
                                name,
                                MoietyType.DEUTERATED_METHYL,
                                moietyType
                            )
                        );
                    }
                }
            });
        }

        /**
         * Sets the list of {@link MFModelIso} model names to evaluate.
         * If not called, all registered models are used (subject to
         * moiety-type filtering).
         *
         * <p>Valid model names are: {@code "1"}, {@code "1f"}, {@code "1s"},
         * {@code "2s"}, {@code "2sf"}, and (only for
         * {@link MoietyType#DEUTERATED_METHYL}) {@code "1sf"}.</p>
         *
         * @param modelNames the model names to evaluate; must not be
         *                   {@code null} or empty, and must contain only
         *                   valid names for the configured moiety type
         * @return this builder
         * @throws NullPointerException if {@code modelNames} is {@code null}
         */
        public Builder modelNames(List<String> modelNames) {
            Objects.requireNonNull(modelNames, "modelNames must not be null");
            this.modelNames = Optional.of(modelNames);
            return this;
        }
    }

    /**
     * Constructs a new {@code ConventionalFitSpec} from the given builder.
     *
     * <p>If no explicit model names were provided, all valid models models are
     * used.</p>
     *
     * @param builder the builder containing all configuration values
     */
    protected ConventionalFitSpec(Builder builder) {
        super(builder);

        List<String> modelNames;
        if (builder.modelNames.isEmpty()) {
            modelNames = new ArrayList<>(Arrays.asList(MFModelIso.getAllModelNames()));
        } else {
            modelNames = new ArrayList<>(builder.modelNames.get());
        }
        switch (moietyType) {
            case AMIDE -> modelNames.remove("1sf");
            case DEUTERATED_METHYL -> {}
        }

        this.modelNames = modelNames;
    }

    /**
     * Serializes this specification to a TOML-formatted string, including
     * the base fields and the list of model names.
     *
     * @return TOML representation of this fit specification
     */
    @Override
    public String toToml() {
        StringBuilder builder = getBaseTomlBuilder();
        builder.append(
            String.format(
                "model_names = [%s]",
                modelNames
                    .stream()
                    .map(s -> String.format("\"%s\"", s))
                    .collect(Collectors.joining(", "))
            )
        );
        return builder.toString();
    }

    /**
     * Returns the list of model names to be evaluated during fitting.
     *
     * @return the model name list (may have been filtered by moiety type)
     */
    public List<String> getModelNames() {
        return modelNames;
    }

    /**
     * Builds all candidate {@link MFModelIso} instances for the given residue.
     *
     * @param data relaxation data for the residue
     * @return list of configured models, one per entry in {@link #getModelNames()}
     * @throws IllegalStateException if tau_M has not been set
     */
    // FIXME: repeated with BaggingFitSpec
    public List<MFModelIso> getModels(MolDataValues<? extends RelaxDataValue> data) {
        return getModelNames()
            .stream()
            .map(name -> getModel(name, data))
            .collect(Collectors.toList());
    }

    /**
     * Appends the model name list to the canonical state string for
     * fingerprinting purposes.
     *
     * @param sb the builder to append to
     */
    @Override
    protected void appendSubclassState(StringBuilder sb) {
        sb.append("modelNames=").append(modelNames.toString()).append('|');
    }

    /**
     * Fits all candidate models to the given residue and selects the best
     * one by AICc.
     *
     * <p>For each model, the method:</p>
     * <ol>
     *   <li>Fits the original (non-resampled) data to obtain the reported
     *       parameter values.</li>
     *   <li>Generates {@code nReplicates} bootstrap samples and fits each
     *       to obtain a series of parameter estimates.</li>
     *   <li>Computes parameter errors from the bootstrap distribution
     *       (via the configured {@link BootstrapMode} strategy).</li>
     *   <li>Constructs an {@link OrderPar} and stores it in
     *       {@code orderParSetMap} under a model-specific key.</li>
     *   <li>Compares the AICc to the current best and updates if improved.</li>
     * </ol>
     *
     * @param key            identifier for the relaxation data set
     * @param data           measured relaxation values for the residue
     * @param orderParSetMap mutable map to which {@link OrderParSet} entries
     *                       for each evaluated model will be added
     * @return the {@link ModelFitResult} for the best model (lowest AICc)
     * @throws AssertionError if no model could be fit (should not happen if
     *                        at least one model is configured)
     * @throws IllegalStateException if tau_M has not been set
     */
    @Override
    public ModelFitResult fit(String key, MolDataValues<?> data, Map<String, OrderParSet> orderParSetMap) {
        RelaxFit relaxFit = initRelaxFit(key, data);
        List<MFModelIso> models = getModels(data);

        Map<String, Score> originalFits = new HashMap<>();
        for (MFModelIso model : models) {
            data.setTestModel(model);
            Score score = runFit(relaxFit, model);
            originalFits.put(model.getName(), score);
        }
        Map.Entry<String, Score> bestOriginalFit = Collections.min(
            originalFits.entrySet(),
            (fit1, fit2) -> Double.compare(fit1.getValue().aicc().get(), fit2.getValue().aicc().get())
        );
        String bestModelName = bestOriginalFit.getKey();
        Score bestScore = bestOriginalFit.getValue();

        int nWeights = data.getNSpectralDensities();
        Optional<ModelFitResult> result =  Optional.empty();
        for (MFModelIso model : models) {
            data.setTestModel(model);
            int nParameters = model.getNPars();

            // Perform bootstrapping to estimate parameter errors
            double[][] parameters = new double[nParameters][nReplicates];
            double[][] weights = new double[nWeights][nReplicates];
            BootstrapSampler<? extends RelaxDataValue> sampler = getBootstrapSampler(data);
            for (int i = 0; i < nReplicates; i++) {
                MolDataValues<? extends RelaxDataValue> replicateData = sampler.sample();
                relaxFit.setRelaxData(key, replicateData);
                Score replicateScore = runFit(relaxFit, model);
                double[] replicateParameters = replicateScore.getPars();
                double[] replicateWeights = replicateData.getWeights();
                for (int k = 0; k < nParameters; k++) parameters[k][i] = replicateParameters[k];
                for (int j = 0; j < nWeights; j++) weights[j][i] = replicateWeights[j];
            }

            String modelName = model.getName();
            double[] fitParameters = originalFits.get(modelName).pars;
            double[] fitErrors = computeStatistics(parameters, weights).getRight();

            String resultKey = makeKey(model.getName());
            orderParSetMap.computeIfAbsent(resultKey, ky -> new OrderParSet(ky));
            makeOrderPar(
                orderParSetMap.get(resultKey),
                sampler.getOriginalData(),
                key,
                originalFits.get(modelName),
                model,
                fitParameters,
                fitErrors
            );

            // Bit hacky, but I'm simply making a duplicate element in the
            // OrderParSetMap for the "optimal" model.
            if (modelName == bestModelName) {
                String bestKey = KEY + "-BEST";
                orderParSetMap.computeIfAbsent(bestKey, ky -> new OrderParSet(ky));
                OrderPar bestOrderPar = makeOrderPar(
                    orderParSetMap.get(bestKey),
                    sampler.getOriginalData(),
                    key,
                    originalFits.get(modelName),
                    model,
                    fitParameters,
                    fitErrors
                );
                result = Optional.of(new ModelFitResult(bestOrderPar, parameters, null));
            }
        }

        return result.orElseThrow(() -> new AssertionError("`result` shouldn't be empty here"));
    }

    /**
     * Constructs the {@link OrderParSet} key for a given model name by
     * combining the strategy prefix with the model identifier.
     *
     * <p>Example: model name {@code "model2sf"} produces key
     * {@code "CONVENTIONAL-2sf"}.</p>
     *
     * @param modelName the model name
     * @return the formatted key string
     */
    private String makeKey(String modelName) {
        return String.format("%s-%s", KEY, modelName.replace("model", ""));
    }
}
