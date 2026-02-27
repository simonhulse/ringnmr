package org.comdnmr.modelfree;

import java.util.Optional;

public class Score {

    final double rss;
    final int nValues;
    final int nPars;
    final boolean parsOK;
    final double complexityS;
    final double complexityTau;
    final double[] pars;

    public Score(double rss, int nValues, int nPars, boolean parsOK, double complexityS, double complexityTau) {
        this(rss, nValues, nPars, parsOK, complexityS, complexityTau, null);
    }

    public Score(double rss, int nValues, int nPars, boolean parsOK, double complexityS, double complexityTau, double[] pars) {
        this.rss = rss;
        this.nValues = nValues;
        this.nPars = nPars;
        this.parsOK = parsOK;
        this.complexityS = complexityS;
        this.complexityTau = complexityTau;
        this.pars = pars;
    }

    public double[] getPars() {
        return pars;
    }

    public double rms() {
        return Math.sqrt(rss / nValues);
    }

    public int getN() {
        return nValues;
    }

    public double value() {
        return value(0.0, 0.0);
    }

    public double value(double lambdaS, double lambdaTau) {
        double score = rms();
        if (!parsOK) {
            score += nValues * 10.0;
        }
        score += complexityS * lambdaS + complexityTau * lambdaTau;
        return score;
    }

    public double complexityS() {
        return complexityS;
    }

    public double complexityTau() {
        return complexityTau;
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
