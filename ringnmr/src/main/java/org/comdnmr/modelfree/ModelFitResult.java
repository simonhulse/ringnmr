package org.comdnmr.modelfree;

import org.nmrfx.chemistry.relax.OrderPar;

public record ModelFitResult(OrderPar orderPar, double[][] replicateData, Double validationValue, Score[] scores) {
}
