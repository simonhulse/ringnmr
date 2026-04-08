package org.comdnmr.modelfree;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.apache.commons.lang3.tuple.Triple;
import org.comdnmr.modelfree.models.MFModelIso;

import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;


public class ConventionalFitSpec extends ModelSelectionFitSpec {

    private static final String KEY = "CONVENTIONAL";

    public static class Builder extends ModelSelectionFitSpec.Builder {

        public ConventionalFitSpec build() {
            return new ConventionalFitSpec(this);
        }

        public Builder modelNames(List<String> modelNames) {
            this.modelNames.clear();
            this.modelNames.addAll(modelNames);
            return this;
        }
    }

    ConventionalFitSpec(Builder builder) { super(builder); }

    @Override
    public ModelFitResult fit(String key, MolDataValues data, Map<String, OrderParSet> orderParSetMap) {
        RelaxFit relaxFit = initRelaxFit(key, data);
        List<MFModelIso> models = getModels(data);

        Optional<Triple<Score, ModelFitResult, MFModelIso>> best = Optional.empty();
        int nWeights = data.getNValues();
        for (MFModelIso model : models) {
            data.setTestModel(model);

            // Perform bootstrapping to estimate parameter errors
            int nParameters = model.getNPars();
            double[][] parameters = new double[nReplicates][nParameters];
            double[][] weights = new double[nReplicates][nWeights];
            BootstrapSampler sampler = getBootstrapSampler(data);

            for (int i = 0; i < nReplicates; i++) {
                MolDataValues replicateData = sampler.sample();
                relaxFit.setRelaxData(key, replicateData);
                Score replicateScore = runFit(relaxFit, model);
                double[] replicateParameters = replicateScore.getPars();
                double[] replicateWeights = replicateData.getWeights();
                for (int k = 0; k < nParameters; k++) parameters[i][k] = replicateParameters[k];
                for (int j = 0; j < nWeights; j++) weights[i][j] = replicateWeights[j];
            }

            // Determine the optimal model using the AICc
            data = sampler.getOriginalData();
            relaxFit.setRelaxData(key, data);
            Score score = runFit(relaxFit, model);

            String resultKey = makeKey(model.getName());
            orderParSetMap.computeIfAbsent(resultKey, k -> new OrderParSet(k));
            OrderPar orderPar = makeOrderParSet(
                orderParSetMap.get(resultKey),
                sampler.getOriginalData(),
                key,
                score,
                model,
                parameters,
                weights
            );

            if (best.isEmpty() || score.aicc().get() < best.get().getLeft().aicc().get()) {
                ModelFitResult result = new ModelFitResult(orderPar, parameters, null);
                best = Optional.of(Triple.of(score, result, model));
            }
        }

        // FIXME: currently don't see the best OrderParSet in RING with this code...
        // String bestModelName = best
        //     .orElseThrow(() -> new AssertionError("`best` should not be empty here!"))
        //     .getRight()
        //     .getName();
        // orderParSetMap.put(makeBestKey(bestModelName), orderParSetMap.get(makeKey(bestModelName)));
        // System.out.printf("orderParSetMap: %s%n", orderParSetMap);

        return best
            .orElseThrow(() -> new AssertionError("`best` should not be empty here!"))
            .getMiddle();
    }

    private String makeKey(String modelName) {
        return String.format("%s-%s", KEY, modelName.replace("model", ""));
    }

    private String makeBestKey(String modelName) {
        return String.format("%s-%s-%s", KEY, "BEST", modelName);
    }
}
