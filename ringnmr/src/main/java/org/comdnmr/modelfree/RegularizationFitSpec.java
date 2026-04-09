package org.comdnmr.modelfree;

import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;

public class RegularizationFitSpec extends FitSpec {

    private static final String KEY = "REGULARIZATION";

    // If s2 is above this number, set to 1 and set corresponding tau to 0
    private static final double S2_THOLD = 0.99;
    // If tau is below 1ps, set to 0s
    private static final double TAU_THOLD = 1.0e-3;

    private final double lambdaS;
    private final double lambdaTauF;
    private final double lambdaTauS;

    public static class Builder extends FitSpec.Builder<Builder> {
        private static final double DEFAULT_LAMBDA_S = 0.5;
        private static final double DEFAULT_LAMBDA_TAUF = 0.1;
        private static final double DEFAULT_LAMBDA_TAUS = 0.2;

        private double lambdaS = DEFAULT_LAMBDA_S;
        private double lambdaTauF = DEFAULT_LAMBDA_TAUF;
        private double lambdaTauS = DEFAULT_LAMBDA_TAUS;

        public static double getDefaultLambdaS() { return DEFAULT_LAMBDA_S; }

        public static double getDefaultLambdaTauF() { return DEFAULT_LAMBDA_TAUF; }

        public static double getDefaultLambdaTauS() { return DEFAULT_LAMBDA_TAUS; }

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

    double getLambdaS() { return lambdaS; }

    double getLambdaTauF() { return lambdaTauF; }

    double getLambdaTauS() { return lambdaTauS; }

    @Override
    public String toToml() {
        StringBuilder builder = getBaseTomlBuilder();
        builder.append(String.format("lambdaS = %s%n", lambdaS));
        builder.append(String.format("lambdaTauF = %s%n", lambdaTauF));
        builder.append(String.format("lambdaTauS = %s", lambdaTauS));
        return builder.toString();
    }

    protected double[] getLower(MFModelIso model) {
        int nParameters = model.getNPars();
        return new double[nParameters];
    }

    protected double[] processParamsAfterFit(MFModelIso model, double[] params) {
        MFModelIso2sf model2sf;
        try {
            model2sf = (MFModelIso2sf) model;
        } catch (ClassCastException exc) {
            throw new AssertionError("Expected regularization model to be 2sf!");
        }
        double s1; double tau1; double s2; double tau2;
        boolean fitTau = model2sf.fitTau();
        int start;
        start = (fitTau) ? 1 : 0;
        s1 = params[start];
        tau1 = params[start + 1];
        s2 = params[start + 2];
        tau2 = params[start + 3];

        double sf2; double tauf; double ss2; double taus;

        // No local motions
        if (s1 > S2_THOLD && s2 > S2_THOLD) {
            // "model 0"
            sf2 = ss2 = 1.0;
            tauf = taus = 0.0;
        }

        // One local motion
        else if (s1 > S2_THOLD || s2 > S2_THOLD) {
            double s;
            double tau;
            if (s1 > S2_THOLD) { s = s2; tau = tau2; }
            else { s = s1; tau = tau1; }
            if (tau < TAU_THOLD) {
                // Model 1
                sf2 = s; tauf = 0.0; ss2 = 1.0; taus = 0.0;
            } else if (tau < model2sf.SLOW_LIMIT) {
                // Model 1f
                sf2 = s; tauf = tau; ss2 = 1.0; taus = 0.0;
            } else {
                // Model 1s
                sf2 = 1.0; tauf = 0.0; ss2 = s; taus = tau;
            }
        }

        // Two local motions
        else {
            if (tau1 < TAU_THOLD) {
                // Model 2s
                sf2 = s1; tauf = 0.0; ss2 = s2; taus = tau2;
            } else if (tau2 < TAU_THOLD) {
                // Model 2s
                sf2 = s2; tauf = 0.0; ss2 = s1; taus = tau1;
            } else {
                // Model 2sf
                if (tau1 < tau2) { sf2 = s1; tauf = tau1; ss2 = s2; taus = tau2; }
                else { sf2 = s2; tauf = tau2; ss2 = s1; taus = tau2; }
            }
        }

        params[start] = sf2;
        params[start + 1] = tauf;
        params[start + 2] = ss2;
        params[start + 3] = taus;

        return params;
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
        double[][] parameters = new double[nParameters][nReplicates];
        double[][] weights = new double[nWeights][nReplicates];
        BootstrapSampler sampler = getBootstrapSampler(data);

        Score[] scores = new Score[nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            Score score = runFit(relaxFit, model);
            scores[i] = score;
            double[] replicateParameters = model.getStandardPars(score.getPars());
            double[] replicateWeights = replicateData.getWeights();
            for (int k = 0; k < nParameters; k++) parameters[k][i] = replicateParameters[k];
            for (int j = 0; j < nWeights; j++) weights[j][i] = replicateWeights[j];
        }

        Pair<double[], double[]> parameterEstimates = computeStatistics(parameters, weights);
        double[] fitParameters = parameterEstimates.getLeft();
        double[] fitErrors = parameterEstimates.getRight();

        orderParSetMap.computeIfAbsent(KEY, ky -> new OrderParSet(ky));
        OrderPar orderPar = makeOrderParSet(
            orderParSetMap.get(KEY),
            sampler.getOriginalData(),
            key,
            // FIXME: not sure what Score should be for makeOrderParSet...
            scores[0],
            model,
            fitParameters,
            fitErrors
        );

        return new ModelFitResult(orderPar, parameters, null);
    }

    protected void appendSubclassState(StringBuilder sb) {
        sb.append("lambdaS=").append(Double.doubleToLongBits(lambdaS)).append('|');
        sb.append("lambdaTauF=").append(Double.doubleToLongBits(lambdaTauF)).append('|');
        sb.append("lambdaTauS=").append(Double.doubleToLongBits(lambdaTauS)).append('|');
    }

}
