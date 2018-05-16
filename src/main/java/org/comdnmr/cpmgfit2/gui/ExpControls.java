/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.gui;

import java.util.List;
import javafx.fxml.FXML;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import org.comdnmr.cpmgfit2.calc.ExpFit;
import org.comdnmr.cpmgfit2.calc.ParValueInterface;
import static org.comdnmr.cpmgfit2.gui.ExpControls.PARS.A;
import static org.comdnmr.cpmgfit2.gui.ExpControls.PARS.C;
import static org.comdnmr.cpmgfit2.gui.ExpControls.PARS.R;

/**
 *
 * @author Bruce Johnson
 */
public class ExpControls implements EquationControls {

    @FXML
    ChoiceBox<String> equationSelector;

    String[] parNames = {"A", "R", "C"};

    enum PARS implements ParControls {
        A("A", 0.0, 500.0, 100.0, 100.0),
        R("R", 0.0, 10.0, 2.0, 2.0),
        C("C", 0.0, 100.0, 20.0, 10.0),;
        String name;
        Slider slider;
        Label label;
        Label valueText;

        PARS(String name, double min, double max, double major, double value) {
            this.name = name;
            slider = new Slider(min, max, value);
            slider.setShowTickLabels(true);
            slider.setShowTickMarks(true);
            slider.setMajorTickUnit(major);
            label = new Label(name);
            label.setPrefWidth(50.0);
            valueText = new Label();
            valueText.setPrefWidth(50);
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
        public void disabled(boolean state) {
            slider.setDisable(state);
        }

        @Override
        public void setValue(double value) {
            slider.setValue(value);
            valueText.setText(String.format("%.1f", value));
        }

        @Override
        public void setText() {
            double value = slider.getValue();
            valueText.setText(String.format("%.1f", value));
        }

        @Override
        public double getValue() {
            return slider.getValue();
        }

    }

    boolean updatingTable = false;
    PyController controller;

    public VBox makeControls(PyController controller) {
        this.controller = controller;
        equationSelector = new ChoiceBox<>();
        equationSelector.getItems().addAll(ExpFit.getEquationNames());
        equationSelector.setValue(ExpFit.getEquationNames().get(0));
        VBox vBox = new VBox();
        HBox hBox1 = new HBox();
        HBox.setHgrow(hBox1, Priority.ALWAYS);
        hBox1.getChildren().add(equationSelector);
        vBox.getChildren().add(hBox1);
        int i = 0;

        for (ParControls control : PARS.values()) {
            HBox hBox = new HBox();
            HBox.setHgrow(hBox, Priority.ALWAYS);
            control.addTo(hBox);

            control.getSlider().valueProperty().addListener(e -> {
                simSliderAction(control.getName());
            });
            vBox.getChildren().add(hBox);
        }

        equationSelector.valueProperty().addListener(e -> {
            equationAction();
        });
        return vBox;
    }

    void equationAction() {
        String equationName = equationSelector.getValue();
        if (equationName == ""){
            equationName = equationSelector.getItems().get(0);
        }
        switch (equationName) {
            case "EXPAB":
                A.disabled(false);
                R.disabled(false);
                C.disabled(true);
                break;
            case "EXPABC":
                A.disabled(false);
                R.disabled(false);
                C.disabled(false);
                break;
            default:
                return;
        }
        // simSliderAction(equationName);

    }

    public void simSliderAction(String label) {
        if (updatingTable) {
            return;
        }
        String equationName = equationSelector.getValue().toString();
        if (equationName.equals("CPMGSLOW") && label.equals("Rex")) {
            return;
        }
        double a = A.getValue();
        double r = R.getValue();
        double c = C.getValue();
        A.setText();
        R.setText();
        C.setText();
        double[] pars;
        switch (equationName) {
            case "EXPAB":
                pars = new double[2];
                pars[0] = a;
                pars[1] = r;

                break;
            case "EXPABC":
                pars = new double[3];
                pars[0] = a;
                pars[1] = r;
                pars[2] = r;

                break;
            default:
                return;
        }
        double[] errs = new double[pars.length];
        int nFields = 1;
        double[] fields = new double[nFields];
        fields[0] = 1.0;
        controller.updateChartEquations(equationName, pars, errs, fields);
    }

    @Override
    public void updateStates(List<int[]> allStates) {

    }

    public void updateSliders(List<ParValueInterface> parValues, String equationName) {
        updatingTable = true;
        for (ParValueInterface parValue : parValues) {
            String parName = parValue.getName();
            ParControls control = PARS.valueOf(parName.toUpperCase());
            if (control != null) {
                control.setValue(parValue.getValue());
            }
        }
        equationSelector.setValue(equationName);
        updatingTable = false;
    }

    public String getEquation() {
        return equationSelector.getValue();
    }

}
