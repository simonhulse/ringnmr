package org.comdnmr.fit;

public record FitQuality(Double rms, Double aic, Double aicc, Double rChiSq, Integer n) {
}
