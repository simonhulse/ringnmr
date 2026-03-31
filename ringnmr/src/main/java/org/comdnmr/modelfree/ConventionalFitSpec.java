package org.comdnmr.modelfree;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.apache.commons.lang3.tuple.Pair;

import org.comdnmr.modelfree.models.MFModelIso;

import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
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

    public ModelFitResult fit(String key, MolDataValues data) {
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setFitJ(fitJ);
        relaxFit.setRelaxData(key, data);

        // Determine the optimal model using the AICc
        List<MFModelIso> models = getModels(data);
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
        data.setTestModel(bestModel);
        int nParameters = bestModel.getNPars();
        int nWeights = data.getNValues();
        double[][] parameters = new double[nReplicates][nParameters];
        double[][] weights = new double[nReplicates][nWeights];
        BootstrapSampler sampler = getBootstrapSampler(data);

        for (int i = 0; i < nReplicates; i++) {
            MolDataValues replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            Score score = runFit(relaxFit, bestModel);
            double[] replicateParameters = score.getPars();
            double[] replicateWeights = replicateData.getWeights();
            for (int k = 0; k < nParameters; k++) parameters[i][k] = replicateParameters[k];
            for (int j = 0; j < nWeights; j++) weights[i][j] = replicateWeights[j];
        }

        MoleculeBase moleculeBase = MoleculeFactory.getActive();
        Map<String, OrderParSet> orderParSetMap = moleculeBase.orderParSetMap();
        orderParSetMap.computeIfAbsent(KEY, ky -> new OrderParSet(ky));
        OrderPar orderPar = makeOrderParSet(
            orderParSetMap.get(KEY),
            sampler.getOriginalData(),
            key,
            bestScore,
            bestModel,
            parameters,
            weights
        );

        ModelFitResult result = new ModelFitResult(
            orderPar,
            parameters,
            null,
            new Score[]{bestScore}
        );

        return result;
    }
}
