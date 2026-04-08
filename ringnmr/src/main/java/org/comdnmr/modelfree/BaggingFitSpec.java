package org.comdnmr.modelfree;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.apache.commons.lang3.tuple.Pair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;

public class BaggingFitSpec extends ModelSelectionFitSpec {

    private static final String KEY = "BAGGING";

    public static class Builder extends ModelSelectionFitSpec.Builder {

        public BaggingFitSpec build() {
            return new BaggingFitSpec(this);
        }
    }

    BaggingFitSpec(Builder builder) { super(builder); }

    @Override
    public ModelFitResult fit(String key, MolDataValues data, Map<String, OrderParSet> orderParSetMap) {
        RelaxFit relaxFit = initRelaxFit(key, data);

        List<MFModelIso> models = getModels(data);
        MFModelIso2sf model2sf = (MFModelIso2sf) getModel("2sf", data);

        int nParameters = 5;
        int nWeights = data.getNValues();
        double[][] parameters = new double[nReplicates][nParameters];
        double[][] weights = new double[nReplicates][nWeights];
        BootstrapSampler sampler = getBootstrapSampler(data);

        Score[] bestScores = new Score[nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            Optional<Pair<MFModelIso, Score>> bestModelScore = Optional.empty();
            for (MFModelIso model : models) {
                data.setTestModel(model);
                Score score = runFit(relaxFit, model);
                if (
                    bestModelScore.isEmpty() ||
                    score.aicc().get() < bestModelScore.get().getRight().aicc().get()
                ) {
                    bestModelScore = Optional.of(Pair.of(model, score));
                }
            }

            if (bestModelScore.isEmpty()) {
                throw new AssertionError("`bestModelScore` should not be empty here!");
            }

            MFModelIso bestModel = bestModelScore.get().getLeft();
            Score bestScore = bestModelScore.get().getRight();
            bestScores[i] = bestScore;
            double[] replicateParameters = bestModel.getStandardPars(bestScore.getPars());
            double[] replicateWeights = replicateData.getWeights();
            for (int k = 0; k < nParameters; k++) parameters[i][k] = replicateParameters[k];
            for (int j = 0; j < nWeights; j++) weights[i][j] = replicateWeights[j];
        }

        orderParSetMap.computeIfAbsent(KEY, ky -> new OrderParSet(ky));
        // FIXME: not sure what Score should be for makeOrderParSet...
        OrderPar orderPar = makeOrderParSet(
            orderParSetMap.get(KEY),
            sampler.getOriginalData(),
            key,
            bestScores[0],
            model2sf,
            parameters,
            weights
        );

        return new ModelFitResult(orderPar, parameters, null);
    }
}
