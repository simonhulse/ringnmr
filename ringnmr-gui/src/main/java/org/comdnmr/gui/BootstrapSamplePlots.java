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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import javafx.beans.Observable;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import org.comdnmr.eqnfit.FitResult;
import org.comdnmr.data.ExperimentResult;
import org.controlsfx.dialog.ExceptionDialog;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.DataSeries;
import org.nmrfx.chart.XYCanvasChart;
import org.nmrfx.chart.XYChartPane;
import org.nmrfx.chart.XYValue;
import org.nmrfx.graphicsio.GraphicsIOException;
import org.nmrfx.graphicsio.SVGGraphicsContext;

/**
 *
 * @author brucejohnson
 */
public class BootstrapSamplePlots {

    PyController pyController;
    Stage stage = null;
    XYCanvasChart activeChart;
    BorderPane borderPane = new BorderPane();
    Scene stageScene = new Scene(borderPane, 500, 500);

    ChoiceBox<String> xArrayChoice = new ChoiceBox<>();
    ChoiceBox<String> yArrayChoice = new ChoiceBox<>();

    public BootstrapSamplePlots(PyController controller) {
        pyController = controller;
    }

    public void showMCplot() {
        //Create new Stage for popup window
        if (stage == null) {
            stage = new Stage();
            Label xlabel = new Label("  X Array:  ");
            Label ylabel = new Label("  Y Array:  ");
            //Populate ChoiceBoxes with fitting variable names
            xArrayChoice.getItems().clear();
            yArrayChoice.getItems().clear();
            try {
                xArrayChoice.valueProperty().addListener((Observable x) -> {
                    updateMCplot();
                });
                yArrayChoice.valueProperty().addListener((Observable y) -> {
                    updateMCplot();
                });
            } catch (NullPointerException npEmc1) {
                Alert alert = new Alert(Alert.AlertType.ERROR);
                alert.setContentText("Error: Fit must first be performed.");
                alert.showAndWait();
                return;
            }
            HBox hBox = new HBox();
            HBox.setHgrow(hBox, Priority.ALWAYS);
            Button exportButton = new Button("Export");
            exportButton.setOnAction(e -> exportBarPlotSVGAction(e));
            Button saveButton = new Button("Save");
            saveButton.setOnAction(e -> saveBootstrapData());

            hBox.getChildren().addAll(exportButton, saveButton, xlabel, xArrayChoice, ylabel, yArrayChoice);
            //Create the Scatter chart
            XYChartPane chartPane = new XYChartPane();
            activeChart = chartPane.getChart();
            borderPane.setTop(hBox);
            borderPane.setCenter(chartPane);
            stage.setScene(stageScene);
        }
        updateMCPlotChoices();
        stage.show();
        updateMCplot();
    }

    void saveBootstrapData() {
        FileWriter writer = null;
        try {
            FileChooser fileChooser = new FileChooser();
            fileChooser.setTitle("Save Bootstrap data");
            fileChooser.setInitialDirectory(pyController.getInitialDirectory());
            File selectedFile = fileChooser.showSaveDialog(null);
            writer = new FileWriter(selectedFile);
            if (selectedFile != null) {
                Map<String, double[]> simsMap = getMCSimsMap();

                int nPar = simsMap.keySet().size();
                double[][] data = new double[nPar][];
                int iPar = 0;
                for (String key : simsMap.keySet()) {
                    if (iPar > 0) {
                        writer.write("\t");
                    }
                    writer.write(key);
                    data[iPar++] = simsMap.get(key);
                }
                for (int i = 0; i < data[0].length; i++) {
                    for (int j = 0; j < nPar; j++) {
                        if (iPar > 0) {
                            writer.write("\t");
                        }
                        writer.write(String.format("%.2f", data[j][i]));
                    }
                    writer.write("\n");
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(BootstrapSamplePlots.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                writer.close();
            } catch (IOException ex) {
                Logger.getLogger(BootstrapSamplePlots.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    void updateMCplot() {
        Map<String, double[]> simsMap = getMCSimsMap();
        if (simsMap != null) {
            Axis xAxis = activeChart.getXAxis();
            Axis yAxis = activeChart.getYAxis();
            String xElem = xArrayChoice.getValue();
            String yElem = yArrayChoice.getValue();
            if ((xElem != null) && (yElem != null)) {
                xAxis.setLabel(xElem);
                yAxis.setLabel(yElem);
                xAxis.setZeroIncluded(false);
                yAxis.setZeroIncluded(false);
                xAxis.setAutoRanging(true);
                yAxis.setAutoRanging(true);
                DataSeries series = new DataSeries();
                activeChart.getData().clear();
                //Prepare XYChart.Series objects by setting data
                series.getData().clear();
                double[] xValues = simsMap.get(xElem);
                double[] yValues = simsMap.get(yElem);
                if ((xValues != null) && (yValues != null)) {
                    for (int i = 0; i < xValues.length; i++) {
                        series.add(new XYValue(xValues[i], yValues[i]));
                    }
                }
                activeChart.getData().add(series);
                activeChart.autoScale(true);
            }
        }
    }

    void updateMCPlotChoices() {
        Map<String, double[]> simsMap = getMCSimsMap();
        xArrayChoice.getItems().clear();
        yArrayChoice.getItems().clear();
        if (simsMap != null) {
            List<String> simParNames = simsMap.keySet().stream().sorted().collect(Collectors.toList());
            xArrayChoice.getItems().addAll(simParNames);
            yArrayChoice.getItems().addAll(simParNames);
            xArrayChoice.setValue(simParNames.get(1));
            yArrayChoice.setValue(simParNames.get(0));
        }
    }

    Map<String, double[]> getMCSimsMap() {
        Map<String, double[]> simsMap = null;
        FitResult showFitResult = null;
        ExperimentResult currentResInfo = pyController.getResidueInfo();
        String currentEquationName = "";
        if (currentResInfo != null) {
            currentEquationName = pyController.getParametersEquation();
            showFitResult = currentResInfo.getFitResult(currentEquationName);
        }
        if (showFitResult == null) {
            currentEquationName = pyController.simControls.getEquation();
            showFitResult = pyController.getFitResult();
        }
        if (showFitResult != null) {
            simsMap = showFitResult.getSimsMap();
        } else {
            System.out.println("no results");
        }
        stage.setTitle("Monte Carlo Results: " + currentEquationName);
        return simsMap;
    }

    @FXML
    void exportBarPlotSVGAction(ActionEvent event) {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Export to SVG");
        fileChooser.setInitialDirectory(pyController.getInitialDirectory());
        File selectedFile = fileChooser.showSaveDialog(null);
        if (selectedFile != null) {
            SVGGraphicsContext svgGC = new SVGGraphicsContext();
            try {
                Canvas canvas = activeChart.getCanvas();
                svgGC.create(canvas.getWidth(), canvas.getHeight(), selectedFile.toString());
                exportChart(svgGC);
                svgGC.saveFile();
            } catch (GraphicsIOException ex) {
                ExceptionDialog eDialog = new ExceptionDialog(ex);
                eDialog.showAndWait();
            }
        }
    }

    protected void exportChart(SVGGraphicsContext svgGC) throws GraphicsIOException {
        svgGC.beginPath();
        activeChart.drawChart(svgGC);
    }

}
