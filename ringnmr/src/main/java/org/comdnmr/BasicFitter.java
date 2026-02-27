package org.comdnmr;

import org.comdnmr.data.ExperimentSet;
import org.comdnmr.eqnfit.ParValueInterface;
import org.nmrfx.chemistry.relax.ResonanceSource;

import java.util.List;

public interface BasicFitter {
    double rms(double[] pars);
    void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues);
    List<ParValueInterface> guessPars(String eqn);
    double[] getSimXDefaults();

    default void setData(ExperimentSet experimentSet, ResonanceSource[] dynSources) {

    }

    default void setupFit(String eqn) {

    }
}
