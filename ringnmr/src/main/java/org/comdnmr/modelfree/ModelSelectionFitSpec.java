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

    public List<String> getModelNames() {
        return modelNames;
    }

    public MFModelIso getModel(String name, MolDataValues data) {
        boolean fitTauM = fitTauM(data);
        return MFModelIso.buildModel(name, fitTauM, tauM, tauFraction, fitExchange);
    }

    public List<MFModelIso> getModels(MolDataValues data) {
        return getModelNames()
            .stream()
            .map(name -> getModel(name, data))
            .collect(Collectors.toList());
    }

}
