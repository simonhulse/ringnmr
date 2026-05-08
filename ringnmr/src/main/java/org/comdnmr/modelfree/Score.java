package org.comdnmr.modelfree;

import java.util.Optional;

public class Score {

    final double rss;
    final int nValues;
    final int nPars;
    final boolean parsOK;
    final double complexityS;
    final double complexityTauF;
    final double complexityTauS;
    final double[] pars;
    protected double[] weights;

    public Score(double rss, int nValues, int nPars, boolean parsOK, double complexityS, double complexityTauF, double complexityTauS) {
        this(rss, nValues, nPars, parsOK, complexityS, complexityTauF, complexityTauS, null);
    }

    public Score(double rss, int nValues, int nPars, boolean parsOK, double complexityS, double complexityTauF, double complexityTauS, double[] pars) {
        this.rss = rss;
        this.nValues = nValues;
        this.nPars = nPars;
        this.parsOK = parsOK;
        this.complexityS = complexityS;
        this.complexityTauF = complexityTauF;
        this.complexityTauS = complexityTauS;
        this.pars = pars;
    }

    public void setWeights(double[] w) { weights = w; }

    public double[] getWeights() { return weights; }

    public double[] getPars() {
        return pars;
    }

    public double rms() {
        double rms = Math.sqrt(rss / nValues);
        return rms;
    }

    public int getN() {
        return nValues;
    }

    public double value() {
        return value(0.0, 0.0, 0.0);
    }

    public double value(double lambdaS, double lambdaTauF, double lambdaTauS) {
        double score = rms();
        if (!parsOK) {
            score += nValues * 10.0;
        }
        score += complexityS * lambdaS + complexityTauF * lambdaTauF + complexityTauS * lambdaTauS;
        return score;
    }

    public double complexityS() {
        return complexityS;
    }

    public double complexityTauF() {
        return complexityTauF;
    }

    public double complexityTauS() {
        return complexityTauS;
    }

    public boolean parsOK() {
        return parsOK;
    }

    public double aic() {
        return 2 * nPars + chiSq();
    }

    public Optional<Double> aicc() {
        return aiccv();
    }

    public Optional<Double> aiccnv() {
        int k = nPars;
        if ((nValues - k -1) < 1) {
            return Optional.empty();
        } else {
            return Optional.of(aic() + 2.0 * k * (k + 1) / (nValues - k - 1));
        }
    }

    public Optional<Double> aiccv() {
        int k = nPars;
        if ((nValues - k) < 1) {
            return Optional.empty();
        } else {
            return Optional.of(aic() + 2.0 * (k + 1) * (k + 2) / (nValues - k));
        }
    }

    public double chiSq() {
        return rss;
    }

    public double reducedChiSq() {
        return rss / (nValues - nPars);
    }
}
