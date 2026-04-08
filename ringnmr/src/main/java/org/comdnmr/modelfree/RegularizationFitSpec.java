package org.comdnmr.modelfree;

import java.util.Map;

import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;

public class RegularizationFitSpec extends FitSpec {

    private static final String KEY = "REGULARIZATION";

    private final double lambdaS;
    private final double lambdaTauF;
    private final double lambdaTauS;

    public static class Builder extends FitSpec.Builder<Builder> {
        private double lambdaS = 0.0;
        private double lambdaTauF = 0.0;
        private double lambdaTauS = 0.0;

        private void validateLambda(String name, double value) {
            if (value < 0.0) {
                throw new IllegalArgumentException(String.format("%s must be >= 0.0", name));
            }
        }

        public Builder lambdaS(double lambdaS) {
            validateLambda("lambdaS", lambdaS);
            this.lambdaS = lambdaS;
            return this;
        }

        public Builder lambdaTauF(double lambdaTauF) {
            validateLambda("lambdaTauF", lambdaTauF);
            this.lambdaTauF = lambdaTauF;
            return this;
        }

        public Builder lambdaTauS(double lambdaTauS) {
            validateLambda("lambdaTauS", lambdaTauS);
            this.lambdaTauS = lambdaTauS;
            return this;
        }

        public RegularizationFitSpec build() {
            return new RegularizationFitSpec(this);
        }
    }

    protected RegularizationFitSpec(Builder builder) {
        super(builder);
        this.lambdaS = builder.lambdaS;
        this.lambdaTauF = builder.lambdaTauF;
        this.lambdaTauS = builder.lambdaTauS;
    }

    @Override
    public String toToml() {
        StringBuilder builder = getBaseTomlBuilder();
        builder.append(String.format("lambdaS = %s%n", lambdaS));
        builder.append(String.format("lambdaTauF = %s%n", lambdaTauF));
        builder.append(String.format("lambdaTauS = %s", lambdaTauS));
        return builder.toString();
    }

    @Override
    protected RelaxFit initRelaxFit(String key, MolDataValues data) {
        RelaxFit relaxFit = super.initRelaxFit(key, data);
        relaxFit.setUseLambda(true);
        relaxFit.setLambdaS(getLambdaS());
        relaxFit.setLambdaTauF(getLambdaTauF());
        relaxFit.setLambdaTauS(getLambdaTauS());
        return relaxFit;
    }

    private MFModelIso2sf getModel2sf(MolDataValues data) {
        boolean fitTauM = fitTauM(data);
        return (MFModelIso2sf) MFModelIso.buildModel("2sf", fitTauM, tauM, tauMFraction, fitExchange);
    }

    @Override
    public ModelFitResult fit(String key, MolDataValues data, Map<String, OrderParSet> orderParSetMap) {
        RelaxFit relaxFit = initRelaxFit(key, data);
        MFModelIso2sf model = getModel2sf(data);
        data.setTestModel(model);

        int nParameters = model.getNPars();
        int nWeights = data.getNValues();
        double[][] parameters = new double[nReplicates][nParameters];
        double[][] weights = new double[nReplicates][nWeights];
        BootstrapSampler sampler = getBootstrapSampler(data);

        Score[] scores = new Score[nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            Score score = runFit(relaxFit, model);
            scores[i] = score;
            double[] replicateParameters = model.getStandardPars(score.getPars());
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
            scores[0],
            model,
            parameters,
            weights
        );

        return new ModelFitResult(orderPar, parameters, null);
    }

    double getLambdaS() { return lambdaS; }

    double getLambdaTauF() { return lambdaTauF; }

    double getLambdaTauS() { return lambdaTauS; }

    @Override
    public int hashCode() {
        int h = 17;
        h = 31 * h + super.hashCode();
        long lambdaSBits = Double.doubleToLongBits(lambdaS);
        h = 31 * h + (int)(lambdaSBits ^ (lambdaSBits >>> 32));
        long lambdaTauFBits = Double.doubleToLongBits(lambdaTauF);
        h = 31 * h + (int)(lambdaTauFBits ^ (lambdaTauFBits >>> 32));
        long lambdaTauSBits = Double.doubleToLongBits(lambdaTauS);
        h = 31 * h + (int)(lambdaTauSBits ^ (lambdaTauSBits >>> 32));
        return h;
    }
}
