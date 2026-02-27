package org.comdnmr.modelfree.models;

public class MFModelSolidState extends MFModelIso1 {

    public MFModelSolidState(boolean fitTau, double targetTau, double tauFraction,
            boolean includeEx) {
        super(fitTau, targetTau, tauFraction, includeEx);
    }

    public MFModelSolidState(double targetTau) { super(targetTau); }

    public MFModelSolidState() { super(); }

    // Spectral density commonly used for solid state (i.e. with no global tumbling term).
    // See Eq 14 of J. Phys. Chem. B 2017, 121, 25, 6117–6130
    @Override
    public double spectralDensity(double s2, double omega, double tau) {
        double oneMinusS2 = 1.0 - s2;
        return super.spectralDensity(oneMinusS2, omega, tau);
    }

    @Override
    public String getName() {
        return "modelSS";
    }
}