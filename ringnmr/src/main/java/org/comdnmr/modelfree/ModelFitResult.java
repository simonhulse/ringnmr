package org.comdnmr.modelfree;

import java.util.Arrays;
import java.util.StringJoiner;

import org.nmrfx.chemistry.relax.OrderPar;

public record ModelFitResult(OrderPar orderPar, double[][] replicateData, Double validationValue, double[] replicateTimes) {

    public String toToml(boolean includeReplicates) {
        StringBuilder builder = new StringBuilder("[fit_result]\n");

        String modelName = orderPar.getModel().replace("model", "").replace("D", "");
        builder.append(String.format("model = \"%s\"%n", modelName));

        // taum
        builder.append(String.format("taum = %s%n", orderPar.getValue("Tau_e")));
        builder.append(String.format("taum_err = %s%n", orderPar.getError("Tau_e")));
        // sf2
        if (modelName.equals("1") || modelName.equals("1f") || modelName.equals("2s") || modelName.equals("2sf")) {
            builder.append(String.format("sf2 = %s%n", orderPar.getValue("Sf2")));
            builder.append(String.format("sf2_err = %s%n", orderPar.getError("Sf2")));
        }
        // tauf
        if (modelName.equals("1f") || modelName.equals("2sf")) {
            builder.append(String.format("tauf = %s%n", orderPar.getValue("Tau_f")));
            builder.append(String.format("tauf_err = %s%n", orderPar.getError("Tau_f")));
        }
        // ss2 & taus
        if (modelName.equals("1s") || modelName.equals("2s") || modelName.equals("2sf")) {
            builder.append(String.format("ss2 = %s%n", orderPar.getValue("Ss2")));
            builder.append(String.format("ss2_err = %s%n", orderPar.getError("Ss2")));
            builder.append(String.format("taus = %s%n", orderPar.getValue("Tau_s")));
            builder.append(String.format("taus_err = %s", orderPar.getError("Tau_s")));
        }

        if (includeReplicates) {
            builder.append("\nreplicates = [\n");
            int nReplicates = replicateData.length;
            for (int i = 0; i < nReplicates; i++) {
                builder.append(String.format("    %s,%n", Arrays.toString(replicateData[i])));
            }
            builder.append("]");
        }

        if (replicateTimes != null) {
            StringJoiner joiner = new StringJoiner(", ", "[", "]");
            for (double t : replicateTimes) {
                joiner.add(String.format("%.6f", t));
            }
            builder.append("\n\n[timing]\n");
            builder.append(String.format("replicate_times_ms = %s", joiner));
        }

        return builder.toString();
    }
}
