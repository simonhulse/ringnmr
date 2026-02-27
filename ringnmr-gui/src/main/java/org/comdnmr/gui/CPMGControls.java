/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.comdnmr.gui;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.control.TextField;
import javafx.scene.input.KeyCode;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import org.comdnmr.data.CPMGExperiment;
import org.comdnmr.data.Experiment;
import org.comdnmr.eqnfit.CPMGFitter;
import org.comdnmr.eqnfit.ParValueInterface;
import org.comdnmr.data.ExperimentResult;
import org.comdnmr.data.ExperimentSet;
import org.comdnmr.eqnfit.CPMGFitFunction;
import org.comdnmr.util.CoMDPreferences;
import org.nmrfx.datasets.Nuclei;

import static org.comdnmr.gui.CPMGControls.PARS.*;

/**
 *
 * @author Bruce Johnson
 */
public class CPMGControls extends EquationControls {

    String[] parNames = {"R2", "Kex", "pA", "delta1", "delta2", "Field2"};

    enum PARS implements ParControls {
        R2("R₂", 0.0, 50.0, 10.0, 10.0, 2),
        KEX("Kₑₓ", 0.0, 20000, 500.0, 500.0, 1),
        PA("pₐ", 0.5, 0.999, 0.1, 0.9, 3),
        // DELTA1 is used for:
        // * Fast regime: \\delta_{ppm}^{min} (see Eq 4 of RING-NMR paper)
        // * Slow regime: \\delta \\omega_H
        // * MQ: \\delta \\omega_H
        DELTA1("", 0.0, 5.0, 0.5, 0.5, 3),
        // DELTA2 is used for:
        // * Fast regime: Nothing
        // * Slow regime: Nothing
        // * MQ: \\delta \\omega_C
        DELTA2("", 0.0, 2.0, 0.5, 0.05, 3),
        FIELD2("B₀", 500.0, 1200.0, 100.0, 600.0, 1),
        TAU("Tau", 0.0, 1.0, 0.2, 0.1, 2)
        ;

        String name;
        Slider slider;
        Label label;
        TextField valueText;
        final String format;

        PARS(String name, double min, double max, double major, double value, int digits) {
            this.name = name;
            slider = new Slider(min, max, value);
            slider.setShowTickLabels(true);
            slider.setShowTickMarks(true);
            slider.setMajorTickUnit(major);
            label = new Label(name);
            label.setPrefWidth(50.0);
            valueText = new TextField();
            valueText.setPrefWidth(70);
            valueText.setText(String.valueOf(value));
            format = "%." + digits + "f";
        }

        @Override
        public String getName() {
            return name;
        }

        @Override
        public void addTo(HBox hBox) {
            hBox.getChildren().addAll(label, slider, valueText);
            HBox.setHgrow(slider, Priority.ALWAYS);
        }

        @Override
        public Slider getSlider() {
            return slider;
        }

        @Override
        public TextField getTextField() {
            return valueText;
        }

        @Override
        public void disabled(boolean state) {
            slider.setDisable(state);
            slider.setVisible(!state);
            valueText.setDisable(state);
            valueText.setVisible(!state);
            label.setVisible(!state);
        }

        @Override
        public void setValue(double value) {
            slider.setValue(value);
            setText();
        }

        @Override
        public void setText() {
            double value = slider.getValue();
            valueText.setText(String.format(format, value));
        }

        @Override
        public double getValue() {
            return slider.getValue();
        }
    }

    boolean updatingTable = false;
    PyController controller;

    @Override
    public VBox makeControls(PyController controller) {
        this.controller = controller;
        VBox vBox = init();
        equationSelector.getItems().addAll(CPMGFitter.getEquationNames());
        equationSelector.setValue(CPMGFitter.getEquationNames().get(0));
        stateSelector.getItems().addAll("0:0:0", "1:0:0");
        stateSelector.setValue("0:0:0");
        HBox hBox1 = new HBox();
        HBox.setHgrow(hBox1, Priority.ALWAYS);
        hBox1.getChildren().addAll(equationSelector, stateSelector, nucleiSelector);
        vBox.getChildren().add(hBox1);

        int i = 0;

        for (ParControls control : PARS.values()) {
            HBox hBox = new HBox();
            HBox.setHgrow(hBox, Priority.ALWAYS);
            control.addTo(hBox);

            control.getSlider().valueProperty().addListener(e -> {
                simSliderAction(control.getName());
                control.setText();
            });

            control.getTextField().setOnKeyReleased(event -> {
                if (event.getCode() == KeyCode.ENTER) {
                    String text = control.getTextField().textProperty().get();
                    if (!text.equals("")) {
                        try {
                            double value = Double.parseDouble(text);
                            control.adjustLimits(value);
                            control.getSlider().setValue(value);
                        } catch (NumberFormatException nfe) {

                        }
                    }
                }
            });

            vBox.getChildren().add(hBox);
        }

        if (controller.simulate == true) {
            equationAction();
        }

        equationSelector.valueProperty().addListener(e -> {
            equationAction();
            updateEquations();
        });
        stateSelector.valueProperty().addListener(e -> {
            stateAction();
        });
        return vBox;
    }

    void equationAction() {
        updatingTable = true;
        String equationName = equationSelector.getValue().toString();
        if (equationName == "") {
            equationName = equationSelector.getItems().get(0);
        }
        ExperimentResult resInfo = controller.chartInfo.getResult();
        if (resInfo != null) {
            updatingTable = true;
            String state = stateSelector.getValue();
            List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
            controller.updateTableWithPars(parValues);
            updateSliders(parValues, equationName);
            ExperimentSet experimentSet = controller.getCurrentExperimentSet();
            if (experimentSet != null) {
                experimentSet.getExperimentData().stream()
                        .findFirst().ifPresent(e -> setNucleus(e.getNucleusName()));
            }
            updatingTable = false;
        }

        switch (equationName) {
            case "NOEX":
                R2.disabled(false);
                KEX.disabled(true);
                PA.disabled(true);
                DELTA1.disabled(true);
                DELTA2.disabled(true);
                FIELD2.disabled(true);
                TAU.disabled(true);
                nucleiSelector.setDisable(false);
                break;

            case "CPMGFAST":
                R2.disabled(false);
                KEX.disabled(false);
                PA.disabled(true);
                DELTA1.disabled(false);
                DELTA1.label.setText("δX");
                DELTA2.disabled(true);
                FIELD2.disabled(true);
                TAU.disabled(true);
                nucleiSelector.setDisable(false);
                break;

            case "CPMGSLOW":
                R2.disabled(false);
                KEX.disabled(false);
                PA.disabled(false);
                DELTA1.disabled(false);
                DELTA1.label.setText("δX");
                DELTA2.disabled(true);
                FIELD2.disabled(true);
                TAU.disabled(true);
                nucleiSelector.setDisable(false);
                break;

            case "CPMGMQ":
                R2.disabled(false);
                KEX.disabled(false);
                PA.disabled(false);
                DELTA1.disabled(false);
                DELTA1.label.setText("δX");
                DELTA2.disabled(false);
                DELTA2.label.setText("δ¹H");
                FIELD2.disabled(true);
                TAU.disabled(false);
                nucleiSelector.setDisable(false);
                break;

            default:
                return;
        }
        simSliderAction(equationName);
        updatingTable = false;
    }

    void stateAction() {
        ExperimentResult resInfo = controller.chartInfo.getResult();
        if (resInfo != null) {
            String state = stateSelector.getValue();
            if (state != null) {
                String equationName = equationSelector.getValue();
                List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
                controller.updateTableWithPars(parValues);
                updateSliders(parValues, equationName);
            }
        }

    }

    // TODO: what is going on here?
    public void simSliderAction(String label) {
        if (updatingTable) {
            return;
        }
        String equationName = equationSelector.getValue().toString();
        if (equationName.equals("CPMGSLOW") && label.equals("dPPMmin")) {
            return;
        }
        updateEquations();
    }

    public double[] getExtras() {
        double B1field = FIELD2.getValue();
        double[] extras = {B1field};
        return extras;
    }

    double[] getPars(String equationName) {
        double r2 = R2.getValue();
        double kEx = KEX.getValue();
        double pA = PA.getValue();
        double delta1 = DELTA1.getValue();
        double delta2 = DELTA2.getValue();
        double field2 = FIELD2.getValue();
        double[] pars;
        switch (equationName) {
            case "NOEX":
                pars = new double[1];
                pars[0] = r2;
                break;
            case "CPMGFAST":
                pars = new double[3];
                pars[0] = kEx;
                pars[1] = r2;
                pars[2] = delta1;  // \\delta_{ppm}^{min}
                break;
            case "CPMGSLOW":
                pars = new double[4];
                pars[0] = kEx;
                pars[1] = pA;
                pars[2] = r2;
                pars[3] = delta1;  // \\delta \\omega_H
                break;
            case "CPMGMQ":
                pars = new double[5];
                pars[0] = kEx;
                pars[1] = pA;
                pars[2] = r2;
                pars[3] = delta1;  // \\delta \\omega_H
                pars[4] = delta2;  // \\delta \\omega_C
                break;
            default:
                pars = null;
        }
        return pars;
    }

    double[] getPars(String equationName, Map<String, ParValueInterface> parValues) {
        double[] pars;
        switch (equationName) {
            case "NOEX":
                pars = new double[1];
                pars[0] = parValues.get("R2").getValue();
                break;
            case "CPMGFAST":
                pars = new double[3];
                pars[0] = parValues.get("Kex").getValue();
                pars[1] = parValues.get("R2").getValue();
                pars[2] = parValues.get("dPPMmin").getValue();
                break;
            case "CPMGSLOW":
                pars = new double[4];
                pars[0] = parValues.get("Kex").getValue();
                pars[1] = parValues.get("pA").getValue();
                pars[2] = parValues.get("R2").getValue();
                pars[3] = parValues.get("dPPM").getValue();
                break;
            case "CPMGMQ":
                pars = new double[5];
                pars[0] = parValues.get("Kex").getValue();
                pars[1] = parValues.get("pA").getValue();
                pars[2] = parValues.get("R2").getValue();
                pars[3] = parValues.get("deltaCPPM").getValue();
                pars[4] = parValues.get("deltaHPPM").getValue();
                break;
            default:
                pars = null;
        }
        return pars;

    }

    public double[] sliderGuess(String equationName, int[][] map) {
        double r2 = R2.getValue();
        double pA = PA.getValue();
        double kEx = KEX.getValue();
        double delta1 = DELTA1.getValue();
        double delta2 = DELTA2.getValue();
        int nPars = CPMGFitFunction.getNPars(map);
        double[] guesses = new double[nPars];
        switch (equationName) {
            case "NOEX":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][0]] = r2;
                }
                break;
            case "CPMGFAST":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][1]] = r2;
                    guesses[map[id][2]] = delta1;  // \\delta_{ppm}^{min}
                }
                guesses[0] = kEx;
                break;
            case "CPMGSLOW":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][2]] = r2;
                    guesses[map[id][3]] = delta1;  // \\delta_X
                }
                guesses[0] = kEx;
                guesses[1] = pA;
                break;
            case "CPMGMQ":
                for (int id = 0; id < map.length; id++) {
                    guesses[map[id][2]] = r2;
                    guesses[map[id][3]] = delta1;  // delta_H
                    guesses[map[id][4]] = delta2;  // delta_C
                }
                guesses[0] = kEx;
                guesses[1] = pA;
                break;
        }
        return guesses;
    }

    @Override
    public void updateStates(List<int[]> allStates) {
        StringBuilder sBuilder = new StringBuilder();
        stateSelector.setDisable(true);
        stateSelector.getItems().clear();
        for (int[] state : allStates) {
            sBuilder.setLength(0);
            for (int i = 1; i < state.length; i++) {
                if (sBuilder.length() > 0) {
                    sBuilder.append(':');
                }
                sBuilder.append(state[i]);
            }
            stateSelector.getItems().add(sBuilder.toString());
        }
        stateSelector.setValue(stateSelector.getItems().get(0));
        stateSelector.setDisable(false);

    }

    @Override
    public void updateSliders(List<ParValueInterface> parValues, String equationName) {
        updatingTable = true;
        for (ParValueInterface parValue : parValues) {
            String parName = parValue.getName();
            if (parName.equalsIgnoreCase("DPPMMIN") || parName.equalsIgnoreCase("DPPM") || parName.equalsIgnoreCase("DELTACPPM")) {
                parName = "DELTA1";
            } else if (parName.equalsIgnoreCase("DELTAHPPM")) {
                parName = "DELTA2";
            }
            ParControls control = PARS.valueOf(parName.toUpperCase());
            if (control != null) {
                control.setValue(parValue.getValue());
            }
        }
        ExperimentSet experimentSet = controller.getCurrentExperimentSet();
        if (experimentSet != null) {
            double[] fields = experimentSet.getFields();
            int iField = Integer.parseInt(stateSelector.getValue().substring(0, 1));
            FIELD2.setValue(fields[iField]);
            FIELD2.setText();
        }
        equationSelector.setValue(equationName);

        updatingTable = false;
    }

    @Override
    public String getEquation() {
        return equationSelector.getValue();
    }

    public List<String> getParNames() {
        String equationName = equationSelector.getValue().toString();
        List<String> parNames1 = new ArrayList<>();
        switch (equationName) {
            case "NOEX":
                parNames1.add("R2");
                parNames1.add("Field2");
                break;
            case "CPMGFAST":
                parNames1.add("Kex");
                parNames1.add("R2");
                parNames1.add("delta1");
                parNames1.add("Field2");
                break;
            case "CPMGSLOW":
                parNames1.add("Kex");
                parNames1.add("pA");
                parNames1.add("R2");
                parNames1.add("delta1");
                parNames1.add("Field2");
                break;
            case "CPMGMQ":
                parNames1.add("Kex");
                parNames1.add("pA");
                parNames1.add("R2");
                parNames1.add("delta1");
                parNames1.add("delta2");
                parNames1.add("Field2");
                break;
        }
        return parNames1;
    }

    void updateEquations() {
        ExperimentResult resInfo = controller.chartInfo.getResult();
        ExperimentSet experimentSet = controller.getCurrentExperimentSet();
        List<GUIPlotEquation> equations = new ArrayList<>();
        double[] pars;
        String equationName = equationSelector.getValue();
        double[] extras = new double[3];
        extras[0] = CoMDPreferences.getRefField() * getNucleus().getFreqRatio();
        extras[1] = CoMDPreferences.getRefField();
        extras[2] = TAU.getValue();
        if (experimentSet == null) {
            pars = getPars(equationName);
            double[] errs = new double[pars.length];
            GUIPlotEquation plotEquation = new GUIPlotEquation("cpmg", equationName, pars, errs, extras);
            equations.add(plotEquation);
        } else {
            String currentState = stateSelector.getValue();
            if (resInfo != null) {
                for (String state : stateSelector.getItems()) {
                    var experimentOpt = experimentSet.getExperimentData().stream().filter(v -> v.getState().equals(state)).findFirst();
                    Experiment experiment = null;
                    if (experimentOpt.isPresent()) {
                        experiment = experimentOpt.get();
                    }
                    Double nucField = null;
                    Double b0Field = null;
                    Double tau = null;
                    if (experiment != null) {
                        nucField = experiment.getNucleusField();
                        b0Field = experiment.getB0Field();
                        if (experiment instanceof CPMGExperiment cpmgExperiment) {
                            tau = cpmgExperiment.getTau();
                        }
                    }
                    int iField = Integer.parseInt(state.substring(0, 1));
                    List<ParValueInterface> parValues = resInfo.getParValues(equationName, state);
                    if (state.equals(currentState) || parValues.isEmpty()) {
                        pars = getPars(equationName);
                    } else {
                        Map<String, ParValueInterface> parMap = new HashMap<>();
                        try {
                            for (ParValueInterface parValue : parValues) {
                                parMap.put(parValue.getName(), parValue);
                            }
                            pars = getPars(equationName, parMap);
                        } catch (NullPointerException npEcpmgpar) {
                            npEcpmgpar.printStackTrace();
                            continue;
                        }
                    }
                    if (nucField != null) {
                        extras[0] = nucField;
                        extras[1] = b0Field;
                        extras[2] = tau;
                    }
                    double[] errs = new double[pars.length];
                    GUIPlotEquation plotEquation = new GUIPlotEquation("cpmg", equationName, pars, errs, extras);
                    equations.add(plotEquation);
                }
            }
        }
        controller.showEquations(equations);
    }
}
