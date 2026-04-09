package org.comdnmr.modelfree;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.comdnmr.modelfree.models.MFModelIso;

public abstract class ModelSelectionFitSpec extends FitSpec {

    protected final List<String> modelNames;

    public static abstract class Builder extends FitSpec.Builder<Builder> {
        protected List<String> modelNames = new ArrayList<>(Arrays.asList("1", "1f", "1s", "2s", "2sf"));
        abstract public ModelSelectionFitSpec build();
    }

    protected ModelSelectionFitSpec(Builder builder) {
        super(builder);
        if (builder.modelNames.isEmpty()) {
            throw new IllegalArgumentException("`modelNames` must not be empty!");
        }
        this.modelNames = builder.modelNames;
    }

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

    public List<String> getModelNames() {
        return modelNames;
    }

    protected double[] getLower(MFModelIso model) { return model.getLower(); }

    protected double[] processParamsAfterFit(MFModelIso model, double[] params) { return params; }

    public MFModelIso getModel(String name, MolDataValues data) {
        boolean fitTauM = fitTauM(data);
        return MFModelIso.buildModel(name, fitTauM, tauM, tauMFraction, fitExchange);
    }

    public List<MFModelIso> getModels(MolDataValues data) {
        return getModelNames()
            .stream()
            .map(name -> getModel(name, data))
            .collect(Collectors.toList());
    }

    protected void appendSubclassState(StringBuilder sb) {
        sb.append("modelNames=").append(modelNames.toString()).append('|');
    }

    @Override
    public int hashCode() {
        int h = 17;
        h = 31 * h + super.hashCode();
        h = 31 * h + (modelNames.hashCode());
        return h;
    }
}
