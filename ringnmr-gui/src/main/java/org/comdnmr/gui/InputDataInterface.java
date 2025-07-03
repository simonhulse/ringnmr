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

import javafx.beans.property.SimpleDoubleProperty;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import org.comdnmr.data.*;
import org.controlsfx.control.textfield.TextFields;
import org.controlsfx.dialog.ExceptionDialog;
import org.nmrfx.chart.DataSeries;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.relax.RelaxTypes;
import org.nmrfx.datasets.DatasetBase;
import org.nmrfx.peaks.PeakList;
import org.nmrfx.utils.GUIUtils;
import org.yaml.snakeyaml.Yaml;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Martha Beckwith
 */
public class InputDataInterface {

    static final String INACTIVE_TEXT_STYLE = "-fx-control-inner-background: red;";
    PyController pyController;

    BorderPane borderPane = new BorderPane();
    GridPane inputInfoDisplay = new GridPane();
    Scene inputScene = new Scene(borderPane);
    Stage infoStage = new Stage();
    TextField chosenDirLabel = new TextField();
    TextField chosenFileLabel = new TextField();
    TextField chosenXPK2FileLabel = new TextField();
    TextField chosenParamFileLabel = TextFields.createClearableTextField();
    ComboBox<String> B0fieldChoice = new ComboBox<>();
    TextField tempTextField = new TextField();
    ChoiceBox<String> nucChoice = new ChoiceBox<>();
    TextField pTextField = new TextField();
    ChoiceBox<String> formatChoice = new ChoiceBox<>();
    TextField tauTextField = new TextField();
    TextArea xValTextArea = new TextArea();
    ChoiceBox<String> fitModeChoice = new ChoiceBox<>();
    ChoiceBox<String> peakListChoice = new ChoiceBox<>();
    TextField B1TextField = new TextField();
    TextField yamlTextField = new TextField();
    CheckBox ppmBox = new CheckBox("ppm to Hz");
    ChoiceBox<String> errModeChoice = new ChoiceBox<>();
    TextField errPercentTextField = new TextField();
    List<HashMap<String, Object>> dataList = new ArrayList<>();
    Button dirChoiceButton = new Button();
    Button fileChoiceButton = new Button();
    Button xpk2ChoiceButton = new Button();
    Button paramFileChoiceButton = new Button();
    Button addButton = new Button();
    Button clearButton = new Button();
    Button yamlButton = new Button();
    Button loadButton = new Button();
    Path dirPath = null;
    ChoiceBox<String> xConvChoice = new ChoiceBox<>();
    ChoiceBox<String> yConvChoice = new ChoiceBox<>();
    TextField delayC0TextField = new TextField();
    TextField delayDeltaTextField = new TextField();
    TextField delayDelta0TextField = new TextField();

    public InputDataInterface(PyController controller) {
        pyController = controller;
    }

    record Choices (ChoiceBox<String> typeChoice,
                    ChoiceBox<DataIO.XCONV> xconvChoiceBox, ChoiceBox<DataIO.YCONV> yconvChoiceBox,
                    TextField tauField, SimpleDoubleProperty tauProperty, ChoiceBox<String> errModeChoiceBox,
                    TextField percentField, SimpleDoubleProperty percentProperty) {

    }

    private void updatePadding() {
        Insets insets = new Insets(15, 15, 15, 15);
        inputInfoDisplay.setPadding(insets);
        inputInfoDisplay.setHgap(10);
        inputInfoDisplay.setVgap(5);
    }
    public void createPeakListInterface() {
        infoStage.setTitle("Load from Peak Lists");
        borderPane.setCenter(inputInfoDisplay);
        inputInfoDisplay.getChildren().clear();
        updatePadding();
        ColumnConstraints col1 = new ColumnConstraints(200);
        ColumnConstraints col2 = new ColumnConstraints(100);
        ColumnConstraints col3 = new ColumnConstraints(150);
        ColumnConstraints col4 = new ColumnConstraints(150);
        ColumnConstraints col5 = new ColumnConstraints(100);
        ColumnConstraints col6 = new ColumnConstraints(120);
        ColumnConstraints col7 = new ColumnConstraints(100);

        double width = 200 + 100 + 150 + 150 + 100 + 120 + 100;
        inputInfoDisplay.setPrefWidth(width);
        borderPane.setPrefWidth(width);
        inputInfoDisplay.getColumnConstraints().clear();
        inputInfoDisplay.getColumnConstraints().addAll(col1, col2, col3, col4, col5, col6, col7);
        List<PeakList> peakLists = new ArrayList<>();
        List<Choices> choices = new ArrayList<>();
        int delta = 1;
        String[] headers = {"PeakList", " Type", "XConv", "YConv", "Tau", "ErrorMode", "Percent"};
        int iCol = 0;
        for (String header : headers) {
            inputInfoDisplay.add(new Label(header), iCol, 0);
            iCol++;
        }
        PeakList.peakLists().forEach(peakList -> {
            if (peakList.hasMeasures()) {
                Label peakListLabel = new Label(peakList.getName());
                peakLists.add(peakList);
                ChoiceBox<String> typeChoice = new ChoiceBox<>();
                typeChoice.getItems().addAll(Arrays.asList("", "R1", "R2", "NOE", "RAP", "RQ", "CPMG"));
                inputInfoDisplay.add(peakListLabel, 0, choices.size() + delta);
                inputInfoDisplay.add(typeChoice, 1, choices.size() + delta);
                ChoiceBox<DataIO.XCONV> xconvChoiceBox = new ChoiceBox<>();
                xconvChoiceBox.getItems().addAll(DataIO.XCONV.values());
                ChoiceBox<DataIO.YCONV> yconvChoiceBox = new ChoiceBox<>();
                yconvChoiceBox.getItems().addAll(DataIO.YCONV.values());
                inputInfoDisplay.add(xconvChoiceBox, 2, choices.size() + delta);
                inputInfoDisplay.add(yconvChoiceBox, 3, choices.size() + delta);
                xconvChoiceBox.setValue(DataIO.XCONV.IDENTITY);
                yconvChoiceBox.setValue(DataIO.YCONV.IDENTITY);
                SimpleDoubleProperty tauProperty = new SimpleDoubleProperty(0.0);
                TextField tauField = GUIUtils.getDoubleTextField(tauProperty);
                ChoiceBox<String> errorModeChoiceBox = new ChoiceBox<>();

                inputInfoDisplay.add(tauField, 4, choices.size() + delta);
                errorModeChoiceBox.getItems().addAll(Arrays.asList("measured", "percent", "replicates"));
                errorModeChoiceBox.setValue("measured");
                inputInfoDisplay.add(errorModeChoiceBox, 5, choices.size() + delta);

                SimpleDoubleProperty percentProperty = new SimpleDoubleProperty(5.0);
                TextField percentField = GUIUtils.getDoubleTextField(percentProperty);
                inputInfoDisplay.add(percentField, 6, choices.size() + delta);

                Choices choices1 = new Choices(typeChoice, xconvChoiceBox, yconvChoiceBox, tauField, tauProperty, errorModeChoiceBox, percentField, percentProperty);
                tauField.setDisable(true);
                percentField.setDisable(true);
                typeChoice.setValue("");
                typeChoice.setOnAction(e -> updateConv(choices1));
                errorModeChoiceBox.setOnAction(e -> updatePercent(choices1));
                choices.add(choices1);

            }
        });
        CheckBox autoFit = new CheckBox("Auto Fit");
        loadButton.setOnAction(e -> loadFromPeakLists(peakLists, choices, autoFit.isSelected()));
        loadButton.setText("Load");
        loadButton.setDisable(false);
        ToolBar toolBar = new ToolBar();
        borderPane.setBottom(toolBar);
        toolBar.getItems().addAll(autoFit, loadButton);
        infoStage.setScene(inputScene);
        infoStage.show();
        infoStage.toFront();
    }

    void updateConv(Choices choices) {
        choices.yconvChoiceBox().setValue(DataIO.YCONV.IDENTITY);
        choices.xconvChoiceBox().setValue(DataIO.XCONV.IDENTITY);
        choices.tauField.setDisable(true);
        switch (choices.typeChoice.getValue()) {
            case "CPMG" -> {
                choices.tauField.setDisable(false);
                choices.yconvChoiceBox().setValue(DataIO.YCONV.RATE);
            }
        }
    }
    void updatePercent(Choices choices) {
        choices.percentField.setDisable(true);
        switch (choices.errModeChoiceBox().getValue()) {
            case "percent" -> {
                choices.percentField.setDisable(false);
            }
        }

    }
    private void loadFromPeakLists(List<PeakList> peakLists, List<Choices> choices, boolean autoFit) {
        for (int i = 0; i < peakLists.size(); i++) {
            Choices choice = choices.get(i);
            String type = choice.typeChoice.getValue();
            if (!type.isBlank()) {
                PeakList peakList = peakLists.get(i);
                ExperimentSet experimentSet = new ExperimentSet(peakList.getName(), peakList.getName());
                if (peakList != null) {
                    int peakDim = 1;
                    String nucleus;
                    double B0field;
                    double temperature;
                    double tau = 0.0;
                    if (type.equalsIgnoreCase("cpmg")) {
                        tau = choice.tauProperty.get();
                        if (tau < 1.0e-6) {
                            GUIUtils.warn("CPMG Tau Value", "Must be non-zero");
                            return;
                        }
                    }
                    DatasetBase dataset = DatasetBase.getDataset(peakList.fileName);
                    if (dataset == null) {
                        nucleus = peakList.getSpectralDim(0).getNucleus();
                        B0field = peakList.getSpectralDim(0).getSf();
                        temperature = 298.14;
                    } else {
                        nucleus = dataset.getNucleus(peakDim).getName();
                        B0field = dataset.getSf(0);
                        temperature = dataset.getTempK();
                    }
                    double percent = choice.percentProperty.doubleValue();
                    DataIO.ErrorMode errorMode = new DataIO.ErrorMode(choice.errModeChoiceBox.getValue(), percent);
                    loadFromPeakList(experimentSet, peakList, type, nucleus, B0field, temperature, tau, choice.xconvChoiceBox.getValue(), choice.yconvChoiceBox.getValue(), errorMode, autoFit);
                }
            }
        }
    }

    public void inputParameters() {

        infoStage.setTitle("Input Data Parameters");
        Label fileLabel = new Label("  Value File:  ");
        Label dirLabel = new Label("  Directory:  ");
        Label peakListLabel = new Label("  PeakList:  ");
        Label xpk2FileLabel = new Label("  XPK2 File:  ");
        Label fitFileLabel = new Label("  CoMD/NMR Analysis File:  ");
        Label fieldLabel = new Label("  B0 Field (1H MHz) :  ");
        Label tempLabel = new Label("  Temperature:  ");
        Label nucLabel = new Label("  Nucleus:  ");
        Label pLabel = new Label("  Pressure:  ");
        Label tauLabel = new Label("  Tau:  ");
        Label xValLabel = new Label("  X Values Conversion:  ");
        Label yValLabel = new Label("  Y Values Conversion:  ");
        Label delayLabel = new Label("  Delays:  ");
        Label fitModeLabel = new Label("  Experiment Type:  ");
        Label B1FieldLabel = new Label("  B1 Field:  ");
        Label yamlLabel = new Label("  YAML File:  ");
        Label errModeLabel = new Label("  Error Mode:  ");
        Label errPercentLabel = new Label("  Error Value:  ");

        Label[] labels = {fitModeLabel, peakListLabel, dirLabel, fileLabel, xpk2FileLabel, fitFileLabel, fieldLabel, tempLabel, pLabel,
                tauLabel, B1FieldLabel, nucLabel, errModeLabel, errPercentLabel, xValLabel, delayLabel, yValLabel, yamlLabel};

        dirChoiceButton.setText("Browse");
        dirChoiceButton.setOnAction(this::chooseDirectory);
        chosenDirLabel.setText("");

        fileChoiceButton.setOnAction(this::chooseFile);
        fileChoiceButton.setText("Browse");
        chosenFileLabel.setText("");
        chosenFileLabel.setStyle("-fx-control-inner-background: red;");
        chosenFileLabel.textProperty().addListener((observable, oldValue, newValue)
                -> {
            if (newValue.equals("")) {
                chosenFileLabel.setStyle("-fx-control-inner-background: red;");
            } else {
                chosenFileLabel.setStyle(null);
            }

        });

        xpk2ChoiceButton.setOnAction(this::chooseXPK2File);
        xpk2ChoiceButton.setText("Browse");
        chosenXPK2FileLabel.setText("");
        chosenXPK2FileLabel.setStyle("-fx-control-inner-background: red;");
        chosenXPK2FileLabel.textProperty().addListener((observable, oldValue, newValue)
                -> {
            if (newValue.equals("")) {
                chosenXPK2FileLabel.setStyle("-fx-control-inner-background: red;");
            } else {
                chosenXPK2FileLabel.setStyle(null);
            }

        });

        paramFileChoiceButton.setOnAction(this::chooseParamFile);
        paramFileChoiceButton.setText("Browse");
        chosenParamFileLabel.setText("");
        chosenParamFileLabel.setStyle("-fx-control-inner-background: red;");
        chosenParamFileLabel.textProperty().addListener((observable, oldValue, newValue)
                -> {
            if (newValue.equals("")) {
                chosenParamFileLabel.setStyle("-fx-control-inner-background: red;");
            } else {
                chosenParamFileLabel.setStyle(null);
            }

        });

        ppmBox.setSelected(false);

        double textFieldWidth = 100;
        double xValAreaWidth = 150; //240;

        TextField[] textFields = {B1TextField, tauTextField, tempTextField, pTextField,
                errPercentTextField};

        for (TextField textField : textFields) {
            textField.setText("");
            textField.setStyle(INACTIVE_TEXT_STYLE);
            textField.textProperty().addListener((observable, oldValue, newValue) -> updateTextField(textField, newValue));
            textField.setMaxWidth(textFieldWidth);
        }

        xValTextArea.setMaxWidth(xValAreaWidth);
        xValTextArea.setWrapText(true);
        yamlTextField.setText("");

        TextField[] texts = {tempTextField, pTextField, tauTextField, B1TextField};

        inputInfoDisplay.getChildren().clear();
        updatePadding();
        borderPane.setCenter(inputInfoDisplay);
        borderPane.setBottom(null);

        for (int i = 0; i < labels.length; i++) {
            inputInfoDisplay.add(labels[i], 0, i);
        }
        for (int i = 0; i < texts.length; i++) {
            inputInfoDisplay.add(texts[i], 1, i + 7);
            texts[i].setMaxWidth(textFieldWidth);
        }

        fitModeChoice.getItems().clear();
        fitModeChoice.getItems().addAll(Arrays.asList("Select", "R1", "R2", "NOE",
                "CPMG", "CEST", "R1RHO", "RAP", "RQ"));
        fitModeChoice.setValue("Select");

        fitModeChoice.valueProperty().addListener(x -> updateInfoInterface());

        peakListChoice.getItems().clear();
        peakListChoice.getItems().add("");
        PeakList.peakLists().forEach(p -> peakListChoice.getItems().add(p.getName()));
        peakListChoice.setValue("");
        peakListChoice.valueProperty().addListener(x -> updatePeakList());

        formatChoice.getItems().clear();
        formatChoice.getItems().addAll(Arrays.asList("mpk2", "ires", "txt"));
        formatChoice.setValue("mpk2");

        nucChoice.getItems().clear();
        nucChoice.getItems().addAll(Arrays.asList("H", "D", "F", "P", "C", "N"));
        nucChoice.setValue("H");

        B0fieldChoice.getItems().clear();
        B0fieldChoice.getItems().addAll(Arrays.asList("400", "475", "500", "600", "700", "750", "800", "900", "950", "1000", "1100", "1200"));
        B0fieldChoice.setValue("");
        B0fieldChoice.valueProperty().addListener((observable, oldValue, newValue)
                -> {
            if (newValue.equals("")) {
                tauTextField.setStyle("-fx-control-inner-background: red;");
            } else {
                tauTextField.setStyle(null);
            }
        });
        B0fieldChoice.setEditable(true);

        EventHandler<ActionEvent> boxevent = e -> {
            String[] xvals = xValTextArea.getText().split("\t");
            ArrayList<Double> fxvals = new ArrayList<>();
            StringBuilder xString = new StringBuilder();
            if ((fitModeChoice.getSelectionModel().getSelectedItem().equals("CEST") || fitModeChoice.getSelectionModel().getSelectedItem().equals("R1RHO"))
                    && ppmBox.isSelected()) {
                for (int i = 0; i < xvals.length; i++) {
                    fxvals.add(Double.parseDouble(xvals[i]) * Double.parseDouble(B0fieldChoice.getSelectionModel().getSelectedItem()));
                    xString.append(fxvals.get(i).toString()).append("\t");
                }
                xValTextArea.setText(xString.toString());
            } else if ((fitModeChoice.getSelectionModel().getSelectedItem().equals("CEST") || fitModeChoice.getSelectionModel().getSelectedItem().equals("R1RHO"))
                    && !ppmBox.isSelected()) {
                for (int i = 0; i < xvals.length; i++) {
                    fxvals.add(Double.parseDouble(xvals[i]) / Double.parseDouble(B0fieldChoice.getSelectionModel().getSelectedItem()));
                    xString.append(fxvals.get(i).toString()).append("\t");
                }
                xValTextArea.setText(xString.toString());
            }
        };

        // set event to checkbox 
        ppmBox.setOnAction(boxevent);

        errModeChoice.getItems().addAll(Arrays.asList("percent", "replicates", "noise", "measured"));
        errModeChoice.setValue("percent");
        errModeChoice.valueProperty().addListener(e -> updateErrorMode());

        xConvChoice.getItems().addAll(Arrays.asList("identity", "tau2", "ppmtohz", "hztoppm", "calc"));
        xConvChoice.setValue("identity");

        xConvChoice.valueProperty().addListener(x -> updateDelays());

        yConvChoice.getItems().addAll(Arrays.asList("identity", "rate", "normalize"));
        yConvChoice.setValue("identity");

        HBox delayBox = new HBox();
        delayBox.getChildren().addAll(new Label("C0:  "), delayC0TextField, new Label("  Delta:  "), delayDeltaTextField, new Label("  Delta0:  "), delayDelta0TextField);

        delayC0TextField.setMaxWidth(textFieldWidth - 20);
        delayDeltaTextField.setMaxWidth(textFieldWidth - 20);
        delayDelta0TextField.setMaxWidth(textFieldWidth - 20);
        int row = 0;
        inputInfoDisplay.add(fitModeChoice, 1, row++);
        inputInfoDisplay.add(peakListChoice, 1, row++);
        inputInfoDisplay.add(dirChoiceButton, 2, row);
        inputInfoDisplay.add(chosenDirLabel, 1, row++);
        inputInfoDisplay.add(fileChoiceButton, 2, row);
        inputInfoDisplay.add(chosenFileLabel, 1, row);
        inputInfoDisplay.add(formatChoice, 3, row++);
        inputInfoDisplay.add(xpk2ChoiceButton, 2, row);
        inputInfoDisplay.add(chosenXPK2FileLabel, 1, row++);
        inputInfoDisplay.add(paramFileChoiceButton, 2, row);
        inputInfoDisplay.add(chosenParamFileLabel, 1, row++);
        inputInfoDisplay.add(B0fieldChoice, 1, row);
        inputInfoDisplay.add(nucChoice, 1, labels.length - 7);
        inputInfoDisplay.add(errModeChoice, 1, labels.length - 6);
        inputInfoDisplay.add(errPercentTextField, 1, labels.length - 5);
        inputInfoDisplay.add(xConvChoice, 1, labels.length - 4);
        inputInfoDisplay.add(delayBox, 1, labels.length - 3, 2, 1);
        inputInfoDisplay.add(yConvChoice, 1, labels.length - 2);
        inputInfoDisplay.add(yamlTextField, 1, labels.length - 1);

        chosenFileLabel.setMaxWidth(200);
        chosenXPK2FileLabel.setMaxWidth(200);
        chosenParamFileLabel.setMaxWidth(200);

        addButton.setOnAction(this::addInfo);
        addButton.setText("Add to Data List");
        inputInfoDisplay.add(addButton, 1, labels.length);

        clearButton.setOnAction(this::clearDataList);
        clearButton.setText("Clear Data List");

        yamlButton.setOnAction(this::makeYAML);
        yamlButton.setText("Create YAML");
        yamlButton.disableProperty().bind(yamlTextField.textProperty().isEmpty());

        loadButton.setOnAction(this::loadInfo);
        loadButton.setText("Load");

        ToolBar toolBar = new ToolBar();
        borderPane.setBottom(toolBar);
        toolBar.getItems().addAll(clearButton, yamlButton, loadButton);

        updateInfoInterface();

        infoStage.setScene(inputScene);
        infoStage.show();
        infoStage.toFront();

    }

    void updateTextField(TextField textField, String newValue) {
        if (newValue.equals("")) {
            textField.setStyle(INACTIVE_TEXT_STYLE);
        } else {
            textField.setStyle(null);
        }
    }

    void updateErrorMode() {
        if ((errModeChoice != null) && (errModeChoice.getValue() != null)) {
            errPercentTextField.setDisable(errModeChoice.getValue().equals("replicates") || errModeChoice.getValue().equals("measured"));
        }
    }

    public void updateInfoInterface() {
        Button[] buttons = {fileChoiceButton, xpk2ChoiceButton, paramFileChoiceButton, addButton, clearButton, loadButton};
        TextField[] textFields = {B1TextField, tauTextField, tempTextField, pTextField,
                errPercentTextField, yamlTextField, chosenFileLabel, chosenXPK2FileLabel, chosenParamFileLabel};
        if (fitModeChoice.getSelectionModel().getSelectedItem() != null) {
            if (fitModeChoice.getSelectionModel().getSelectedItem().equals("Select")) {
                for (TextField textField : textFields) {
                    textField.setDisable(true);
                }
                for (Button button : buttons) {
                    button.setDisable(true);
                }
                ppmBox.setDisable(true);
                xValTextArea.setDisable(true);
                errModeChoice.setDisable(true);
                xConvChoice.setDisable(true);
                yConvChoice.setDisable(true);
                formatChoice.setDisable(true);
                nucChoice.setDisable(true);
                B0fieldChoice.setDisable(true);
                delayC0TextField.setDisable(true);
                delayDeltaTextField.setDisable(true);
                delayDelta0TextField.setDisable(true);
            } else if (!fitModeChoice.getSelectionModel().getSelectedItem().equals("Select")) {
                for (TextField textField : textFields) {
                    textField.setDisable(false);
                }
                for (Button button : buttons) {
                    button.setDisable(false);
                }
                formatChoice.setDisable(false);

                if (!peakListChoice.getValue().equals("")) {
                    TextField[] textFields2 = {yamlTextField, chosenFileLabel, chosenXPK2FileLabel, chosenParamFileLabel};
                    for (TextField textField : textFields2) {
                        textField.setDisable(true);
                    }
                    for (Button button : buttons) {
                        button.setDisable(true);
                    }
                    dirChoiceButton.setDisable(true);
                    formatChoice.setDisable(true);
                    loadButton.setDisable(false);
                }
                errModeChoice.getItems().setAll(Arrays.asList("percent", "replicates", "noise", "measured"));
                errModeChoice.setValue("percent");
                xValTextArea.setDisable(false);
                errModeChoice.setDisable(false);
                xConvChoice.setDisable(false);
                yConvChoice.setDisable(false);
                formatChoice.setDisable(false);
                nucChoice.setDisable(false);
                B0fieldChoice.setDisable(false);
                if (fitModeChoice.getSelectionModel().getSelectedItem().equals("CPMG")) {
                    B1TextField.setDisable(true);
                    tauTextField.setDisable(false);
                    ppmBox.setDisable(true);
                    xConvChoice.getItems().clear();
                    xConvChoice.getItems().addAll(Arrays.asList("identity", "tau2"));
                    yConvChoice.getItems().clear();
                    yConvChoice.getItems().addAll(Arrays.asList("identity", "rate"));
                    xConvChoice.setValue("identity");
                    yConvChoice.setValue("rate");
                } else if (fitModeChoice.getSelectionModel().getSelectedItem().equals("EXP")) {
                    B1TextField.setDisable(true);
                    tauTextField.setDisable(true);
                    ppmBox.setDisable(true);
                    xConvChoice.getItems().clear();
                    xConvChoice.getItems().addAll(Arrays.asList("identity", "calc"));
                    yConvChoice.getItems().clear();
                    yConvChoice.getItems().addAll(Arrays.asList("identity", "normalize"));
                    xConvChoice.setValue("identity");
                    yConvChoice.setValue("identity");
                } else if (fitModeChoice.getSelectionModel().getSelectedItem().equals("NOE")) {
                    errModeChoice.getItems().setAll(Arrays.asList("percent", "noise", "measured"));
                    errModeChoice.setValue("noise");
                    B1TextField.setDisable(true);
                    tauTextField.setDisable(true);
                    ppmBox.setDisable(true);
                    xConvChoice.getItems().clear();
                    xConvChoice.getItems().addAll(List.of("identity"));
                    yConvChoice.getItems().clear();
                    yConvChoice.getItems().addAll(List.of("normalize"));
                    xConvChoice.setValue("identity");
                    yConvChoice.setValue("normalize");
                } else if ((fitModeChoice.getSelectionModel().getSelectedItem().equals("CEST") || fitModeChoice.getSelectionModel().getSelectedItem().equals("R1RHO"))) {
                    B1TextField.setDisable(false);
                    tauTextField.setDisable(false);
                    ppmBox.setDisable(false);
                    xConvChoice.getItems().clear();
                    xConvChoice.getItems().addAll(Arrays.asList("identity", "ppmtohz", "hztoppm"));
                    yConvChoice.getItems().clear();
                    yConvChoice.getItems().addAll(Arrays.asList("identity", "normalize"));
                    xConvChoice.setValue("identity");
                    yConvChoice.setValue("identity");
                    if (fitModeChoice.getSelectionModel().getSelectedItem().equals("CEST")) {
                        yConvChoice.setValue("normalize");
                    }
                }
            }
        }

    }

    public void updateDelays() {
        if (!fitModeChoice.getSelectionModel().getSelectedItem().equals("Select")
                && (xConvChoice.getSelectionModel().getSelectedItem() != null) && xConvChoice.getSelectionModel().getSelectedItem().equals("calc")) {
            delayC0TextField.setDisable(false);
            delayDeltaTextField.setDisable(false);
            delayDelta0TextField.setDisable(false);
        } else {
            delayC0TextField.setDisable(true);
            delayDeltaTextField.setDisable(true);
            delayDelta0TextField.setDisable(true);
        }
    }

    public void chooseDirectory(ActionEvent event) {
        DirectoryChooser dirChooser = new DirectoryChooser();
        File file = dirChooser.showDialog(infoStage);
        if (file != null) {
            chosenDirLabel.setText(file.toString());
            dirPath = file.toPath();
        }
    }

    public void chooseFile(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        if (dirPath != null) {
            File userDirectory = new File(dirPath.toString());
            fileChooser.setInitialDirectory(userDirectory);
        }
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("mpk2 or txt File", "*.mpk2", "*.txt"));
        File file = fileChooser.showOpenDialog(infoStage);
        if (file != null) {
            Path path = dirPath.relativize(file.toPath());
            String pathString = path.toString();
            chosenFileLabel.setText(pathString);
            if (pathString.contains(".")) {
                String fileTail = pathString.substring(0, pathString.lastIndexOf('.'));
                yamlTextField.setText(fileTail + ".yaml");
                if (pathString.endsWith(".mpk2")) {
                    FileSystem fileSystem = FileSystems.getDefault();
                    File xpk2File = fileSystem.getPath(file.getParent(), fileTail + ".xpk2").toFile();
                    if (xpk2File.canRead()) {
                        parseXPK2File(xpk2File);
                    }
                }
            }
        }
    }

    public void chooseXPK2File(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        if (dirPath != null) {
            File userDirectory = new File(dirPath.toString());
            fileChooser.setInitialDirectory(userDirectory);
        }
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("xpk2 File", "*.xpk2"));
        File file = fileChooser.showOpenDialog(infoStage);
        if (file != null) {
            parseXPK2File(file);
        }
    }

    void parseXPK2File(File file) {
        Path path = dirPath.relativize(file.toPath());
        chosenXPK2FileLabel.setText(path.toString());

        Path path1 = file.toPath();

        List<String[]> head = new ArrayList<>();

        try (BufferedReader fileReader = Files.newBufferedReader(path1)) {
            while (true) {
                String line = fileReader.readLine();
                if (line == null) {
                    break;
                }
                String sline = line.trim();
                if (sline.length() == 0) {
                    continue;
                }
                if (sline.startsWith("id")) {
                    break;
                }
                String[] sline1 = line.split("\t", -1);
                head.add(sline1);
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        int sfInd = Arrays.asList(head.get(2)).indexOf("sf");
        int codeInd = Arrays.asList(head.get(2)).indexOf("code");
        String field = Arrays.asList(head.get(3)).get(sfInd);
        String nuc = Arrays.asList(head.get(4)).get(codeInd);
        String nuc1 = nuc.replaceAll("[^a-zA-Z]", "");
        nucChoice.setValue(nuc1);
        B0fieldChoice.getSelectionModel().select(field);
    }

    void updatePeakList() {
        updateInfoInterface();
        PeakList peakList = PeakList.get(peakListChoice.getValue());
        if (peakList != null) {
            int peakDim = 1;
            String nucleus;
            double B0field;
            double temperature;
            DatasetBase dataset = DatasetBase.getDataset(peakList.fileName);
            if (dataset == null) {
                nucleus = peakList.getSpectralDim(0).getNucleus();
                B0field = peakList.getSpectralDim(0).getSf();
                temperature = 298.14;
            } else {
                nucleus = dataset.getNucleus(peakDim).getName();
                B0field = dataset.getSf(0);
                temperature = dataset.getTempK();
                System.out.println(temperature);
            }

            nucChoice.setValue(nucleus);
            B0fieldChoice.setValue(String.valueOf(B0field));
            tempTextField.setText(String.valueOf(temperature));

        }
    }

    Double getDouble(String str) {
        Double value;
        try {
            value = Double.parseDouble(str);
        } catch (NumberFormatException nfe) {
            value = null;
        }
        return value;
    }

    void loadFromPeakList() {
        PeakList peakList = PeakList.get(peakListChoice.getValue());
        String expMode = fitModeChoice.getValue();
        double b0Field = Double.parseDouble(B0fieldChoice.getValue());
        Double temperatureK = getDouble(tempTextField.getText());
        Double tau = getDouble("tau");
        String nucName = nucChoice.getValue();
        DataIO.XCONV xConv = xConvChoice.getValue() == null ? DataIO.XCONV.IDENTITY : DataIO.XCONV.valueOf(xConvChoice.getValue());
        DataIO.YCONV yConv = yConvChoice.getValue() == null ? DataIO.YCONV.IDENTITY :  DataIO.YCONV.valueOf(yConvChoice.getValue());
        if (expMode.equalsIgnoreCase("noe")) {
            yConv = DataIO.YCONV.NORMALIZE;
        } else if (expMode.equalsIgnoreCase("cpmg")) {
            yConv = DataIO.YCONV.RATE;
            xConv = DataIO.XCONV.TAU4;
        }
        DataIO.ErrorMode errorMode = new DataIO.ErrorMode(errModeChoice.getValue(), Double.parseDouble(errPercentTextField.getText()));
        ExperimentSet experimentSet = new ExperimentSet("peaks", "peaks");
        loadFromPeakList(experimentSet, peakList, expMode, nucName, b0Field, temperatureK, tau, xConv, yConv, errorMode, false);
    }

    void loadFromPeakList(ExperimentSet experimentSet, PeakList peakList, String expMode, String nucName, double b0Field, double temperatureK, double tau,
                          DataIO.XCONV xConv, DataIO.YCONV yConv, DataIO.ErrorMode errorMode,
                          boolean autoFit) {
        expMode = expMode.toLowerCase();
        if (peakList != null) {

            MoleculeBase mol = MoleculeFactory.getActive();
            boolean dynCreateMol;
            if (mol == null) {
                if (!GUIUtils.affirm("No molecule present, dynamically create")) {
                    return;
                } else {
                    dynCreateMol = true;
                }
            } else {
                dynCreateMol = false;
            }
            boolean dynCreateAtom = true;
            DynamicsSource dynSource = new DynamicsSource(false, false, dynCreateMol, dynCreateAtom);
            String peakListName = peakList.getName();
            experimentSet.setExpMode(expMode);

            Experiment expData;
            switch (expMode) {
                case "rq":
                case "rap":
                case "r1":
                    expData = new T1Experiment(experimentSet, peakList.getName(),
                            nucName, b0Field, temperatureK);
                    break;
                case "r2":
                    expData = new T2Experiment(experimentSet, peakList.getName(),
                            nucName, b0Field, temperatureK);
                    break;
                case "noe":
                    expData = new NOEExperiment(experimentSet, peakList.getName(),
                            nucName, b0Field, temperatureK);
                case "cpmg":
                    expData = new CPMGExperiment(experimentSet, peakList.getName(),
                            nucName, b0Field, tau, temperatureK);
                    break;
                default:
                    expData = new Experiment(experimentSet, peakList.getName(),
                            nucName, b0Field, temperatureK, expMode);
            }

            try {
                DataIO.loadFromPeakList(peakList, expData, experimentSet,
                        xConv, yConv, errorMode, dynSource);
                ResidueChart reschartNode = PyController.mainController.getActiveChart();
                if (reschartNode == null) {
                    reschartNode = PyController.mainController.addChart();

                }
                reschartNode.getData().clear();
                ChartUtil.addResidueProperty(experimentSet.name(), experimentSet);
                String parName = "Kex";
                if (expMode.equalsIgnoreCase("r1") || expMode.equalsIgnoreCase("r2") ||
                        expMode.equalsIgnoreCase("rq") || expMode.equalsIgnoreCase("rap")) {
                    parName = "R";
                } else if (experimentSet.getExpMode().equalsIgnoreCase("noe")) {
                    parName = "NOE";
                }
                ObservableList<DataSeries> data = ChartUtil.getParMapData(experimentSet.name(), "best", "0:0:0", parName);
                PyController.mainController.setCurrentExperimentSet(experimentSet);
                PyController.mainController.makeAxisMenu();
                PyController.mainController.setYAxisType(experimentSet.getExpMode(), experimentSet.name(),
                        "best", "0:0:0", parName, true);
                reschartNode.setResProps(experimentSet);
                PyController.mainController.setControls();
                if (expMode.equalsIgnoreCase("noe")) {
                    DataIO.addRelaxationFitResults(experimentSet, RelaxTypes.NOE);
                } else {
                    if (autoFit) {
                        PyController.mainController.fitResiduesNow();
                    }
                }
            } catch (IllegalArgumentException iAE) {
                GUIUtils.warn("Load from peak list", iAE.getMessage());
            }
        }
    }

    public void chooseParamFile(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        if (dirPath != null) {
            File userDirectory = new File(dirPath.toString());
            fileChooser.setInitialDirectory(userDirectory);
        }
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("txt File", "*.txt"));
        File file = fileChooser.showOpenDialog(infoStage);
        if (file != null) {
            String directory = file.getParent();
            String fileName = file.getName();
            chosenParamFileLabel.setText(directory + "/" + fileName);
        }
    }

    public void clearParamFile(ActionEvent event) {
        chosenParamFileLabel.setText("");
    }

    public void addInfo(ActionEvent event) {
        addInfo();
    }

    public void addInfo() {
        HashMap<String, Object> hm = new HashMap<>();
        hm.put("file", chosenFileLabel.getText());

        hm.put("paramFile", chosenParamFileLabel.getText());
        hm.put("temperature", Double.parseDouble(tempTextField.getText()));
        hm.put("xconv", xConvChoice.getSelectionModel().getSelectedItem());
        hm.put("yconv", yConvChoice.getSelectionModel().getSelectedItem());
        hm.put("B0", Double.parseDouble(B0fieldChoice.getSelectionModel().getSelectedItem()));
        hm.put("nucleus", nucChoice.getSelectionModel().getSelectedItem().replaceAll("[^a-zA-Z]", ""));
        if (!tauTextField.isDisabled()) {
            hm.put("tau", Double.parseDouble(tauTextField.getText()));
        }
        hm.put("pressure", Double.parseDouble(pTextField.getText()));
        hm.put("format", formatChoice.getSelectionModel().getSelectedItem());
        hm.put("fitmode", fitModeChoice.getSelectionModel().getSelectedItem().toLowerCase());
        if (!B1TextField.isDisabled()) {
            hm.put("B1", Double.parseDouble(B1TextField.getText()));
        }
        HashMap<String, java.io.Serializable> hmde = new HashMap<>();
        hmde.put("mode", errModeChoice.getValue());
        if (!errModeChoice.getValue().equals("replicates") && !errModeChoice.getValue().equals("measured")) {
            hmde.put("value", Double.parseDouble(errPercentTextField.getText()));
        }
        HashMap<String, Double> hmdd = new HashMap<String, Double>();
        if (!delayC0TextField.getText().equals("") && !delayDeltaTextField.getText().equals("") && !delayDelta0TextField.getText().equals("")) {
            hmdd.put("c0", Double.parseDouble(delayC0TextField.getText()));
            hmdd.put("delta0", Double.parseDouble(delayDelta0TextField.getText()));
            hmdd.put("delta", Double.parseDouble(delayDeltaTextField.getText()));
        }

        hm.put("error", hmde);
        if (hm.get("xconv").equals("calc")) {
            hm.put("delays", hmdd);
        }

        String[] xvals = xValTextArea.getText().trim().split("\t");
        if (xvals.length > 0) {
            ArrayList<Double> fxvals = new ArrayList<>();
            try {
                for (String xval : xvals) {
                    fxvals.add(Double.parseDouble(xval));
                }
            } catch (NumberFormatException nfe) {
                fxvals = null;
            }
            hm.put("vcpmg", fxvals);
        }
        dataList.add(hm);
    }

    public void makeYAML(ActionEvent event) {
        if (dataList.isEmpty()) {
            addInfo();
        }
        makeYAML(dataList);
        dataList.clear();
    }

    public void makeYAML(List data) {
        HashMap hm1 = new HashMap();
        HashMap hm2 = new HashMap();
        List<HashMap<String, Object>> dataHmList = new ArrayList<>();

        for (Object datum : data) {
            HashMap hmdf = (HashMap) datum;
            HashMap hmd = new HashMap(hmdf);

            String paramFile = (String) hmdf.get("paramFile");
            String paramFileName = paramFile.substring(paramFile.lastIndexOf("/") + 1);
            hm2.put("mode", hmdf.get("fitmode"));
            if (!paramFileName.equals("")) {
                hm2.put("file", paramFileName);
            } else {
                hm2.put("file", "analysis.txt");
            }
            Set keySet = hmd.keySet();
            if (!hmd.get("fitmode").equals("cest") && !hmd.get("fitmode").equals("r1rho")) {
                keySet.remove("B1");
            }
            if ((hmd.get("vcpmg") == null) || (hmd.get("vcpmg").toString().equals(""))) {
                keySet.remove("vcpmg");
            }
            if (!hmd.get("xconv").equals("calc")) {
                keySet.remove("delays");
            }
            keySet.remove("fitmode");
            keySet.remove("paramFile");
            hmd.keySet().retainAll(keySet);
            String dataFile = (String) hmdf.get("file");
            String dataFileName = dataFile.substring(dataFile.lastIndexOf("/") + 1, dataFile.length());
            hmd.put("file", dataFileName);
            dataHmList.add(hmd);
        }
        hm2.put("data", dataHmList);
        hm1.put("fit", hm2);

        Yaml yaml = new Yaml();
        String s = yaml.dumpAsMap(hm1);
        Path path = FileSystems.getDefault().getPath(dirPath.toString(), yamlTextField.getText());
        try (FileWriter writer = new FileWriter(path.toFile())) {
            writer.write(s);
            System.out.println(yamlTextField.getText() + " written");
        } catch (IOException ex) {
            Logger.getLogger(DataIO.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void clearDataList(ActionEvent event) {
        dataList.clear();
    }

    public void loadInfo(ActionEvent event) {
        if (!peakListChoice.getValue().equals("")) {
            loadFromPeakList();
            return;
        }
        if (dataList.isEmpty()) {
            addInfo();
        }
        String projectName = yamlTextField.getText().trim();
        if (projectName.length() == 0) {

        } else if (projectName.endsWith(".yaml")) {
            projectName = projectName.substring(0, projectName.indexOf(".yaml"));
        }
        File projectDirFile = new File(chosenDirLabel.getText().trim());
        dirPath = projectDirFile.toPath();

        ExperimentSet resProp;
        String expMode = fitModeChoice.getSelectionModel().getSelectedItem().toLowerCase();
        resProp = new ExperimentSet(projectName, projectDirFile.toString());
        expMode = expMode.toLowerCase();
        resProp.setExpMode(expMode);

        try {
            DataIO.processYAMLDataSections(resProp, dirPath, expMode, dataList);
        } catch (IOException ex) {
            ExceptionDialog dialog = new ExceptionDialog(ex);
            dialog.showAndWait();
            return;
        }

        PyController.mainController.clearSecondaryStructure();
        if (PyController.mainController.activeChart != null) {
            PyController.mainController.clearChart();
            PyController.mainController.chartInfo.clear();
            PyController.mainController.simulate = false;
            PyController.mainController.fitResult = null;
        }
        ResidueChart reschartNode = PyController.mainController.getActiveChart();
        if (reschartNode == null) {
            reschartNode = PyController.mainController.addChart();

        }
//            ExperimentSet resProp = DataIO.loadParameters(fileName);
        ChartUtil.addResidueProperty(resProp.name(), resProp);
        String parName = "Kex";
        if (expMode.equals("r1") || expMode.equals("r2") || expMode.equals("rq") || expMode.equals("rap")) {
            parName = "R";
        }
        ObservableList<DataSeries> data = ChartUtil.getParMapData(resProp.name(), "best", "0:0:0", parName);
        PyController.mainController.setCurrentExperimentSet(resProp);
        PyController.mainController.makeAxisMenu();
        PyController.mainController.setYAxisType(resProp.getExpMode(), resProp.name(),
                "best", "0:0:0", parName, true);
        reschartNode.setResProps(resProp);
        PyController.mainController.setControls();

        fitModeChoice.setValue("Select");
        dataList.clear();
        updateInfoInterface();
    }

}
