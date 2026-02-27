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

import java.io.BufferedReader;

import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.canvas.Canvas;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;

import javafx.collections.ListChangeListener;
import javafx.geometry.Orientation;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.input.MouseEvent;

import javafx.scene.paint.Color;
import javafx.stage.FileChooser;
import org.comdnmr.eqnfit.PlotEquation;
import org.comdnmr.data.ExperimentData;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.DataSeries;

import org.nmrfx.chart.XYCanvasChart;
import org.nmrfx.chart.XYValue;
import org.nmrfx.graphicsio.GraphicsContextInterface;
import org.nmrfx.graphicsio.GraphicsContextProxy;
import org.nmrfx.graphicsio.GraphicsIOException;
import org.nmrfx.graphicsio.SVGGraphicsContext;

public class PlotData extends XYCanvasChart {

    ObservableList<DataSeries> simData = FXCollections.observableArrayList();
    String fileName;
    ObservableList<GUIPlotEquation> plotEquations = FXCollections.observableArrayList();

    protected static final Color[] colors = {
            Color.web("#1b9e77"),
            Color.web("#d95f02"),
            Color.web("#7570b3"),
            Color.web("#e7298a"),
            Color.web("#66a61e"),
            Color.web("#e6ab02"),
            Color.web("#a6761d"),
            Color.web("#666666"),
            Color.web("#ff7f00"),
            Color.web("#6a3d9a"),};

    public PlotData(Canvas canvas, final Axis... axes) {
        super(canvas, axes);
        init();
    }

    public static PlotData buildChart(Canvas canvas) {
        Axis xAxis = new Axis(Orientation.HORIZONTAL, 0, 100, 400, 100.0);
        Axis yAxis = new Axis(Orientation.VERTICAL, 0, 100, 100, 400);
        return new PlotData(canvas, xAxis, yAxis);
    }

    void init() {
        xAxis.setUpperBound(1000.0);
        yAxis.setUpperBound(60.0);
        yAxis.setZeroIncluded(true);
        xAxis.setLabel("Time (s)");
        yAxis.setLabel("Intensity");
        xAxis.setLabel("ν (cpmg)");
        yAxis.setLabel("R2 (ν)");
        plotEquations.addListener((ListChangeListener) (e -> drawChart()));
        getCanvas().setOnMouseClicked(this::mouseClicked);
    }

    @Override
    public void setNames(String title, String xName, String yName, String yPad) {
        setTitle(title);
        xAxis.setLabel(xName);
        yAxis.setLabel(yName);
    }

    @Override
    public void setBounds(double xLower, double xUpper, double yLower, double yUpper, double xtick, double ytick) {
        xAxis.setLowerBound(xLower);
        xAxis.setUpperBound(xUpper);
        yAxis.setLowerBound(yLower);
        yAxis.setUpperBound(yUpper);
        drawChart();
    }

    void mouseClicked(MouseEvent e) {
        Optional<Hit> hitOpt = pickChart(e.getX(), e.getY(), 5);
        if (hitOpt.isPresent()) {
            Hit hit = hitOpt.get();
            PyController.mainController.selectTableRow(hit.getSeries().getName(), hit.getIndex());
            PyController.mainController.statusBar.setText(hit.toString());

        }
    }

    public int getNumPlots() {
        return getData().size();
    }

    public int getNumEquations() {
        return plotEquations.size();
    }

    public void setEquations(List<GUIPlotEquation> plotEquations) {
        this.plotEquations.clear();
        this.plotEquations.addAll(plotEquations);
        drawChart();
    }

    public void clear() {
        plotEquations.clear();
        getData().clear();
    }

    @Override
    public void drawChart() {
        super.drawChart();
        paintLines();
    }

    protected void exportVectorGraphics(SVGGraphicsContext svgGC) throws GraphicsIOException {
        svgGC.save();
        svgGC.clip();
        svgGC.beginPath();
        drawChart(svgGC);
        paintLines(svgGC);
        svgGC.restore();
    }

    void dumpNode(String seriesName, XYValue value) {
        Object extraValue = value.getExtraValue();
        if (extraValue instanceof ExperimentData.DataValue dataValue) {
            PyController.mainController.selectTableRow(seriesName, dataValue.getIndex());
        }
    }

    public double[] getXBounds() {
        return new double[]{xAxis.getLowerBound(), xAxis.getUpperBound()};
    }

    public double[] getYBounds() {
        return new double[]{yAxis.getLowerBound(), yAxis.getUpperBound()};
    }

    void paintLines() {
        GraphicsContext gCC = getCanvas().getGraphicsContext2D();
        GraphicsContextInterface gC = new GraphicsContextProxy(gCC);
        paintLines(gC);
    }

    void paintLines(GraphicsContextInterface gC) {

        int nIncr = 256;
        double[] xValues = new double[nIncr];
        double[] yValues = new double[nIncr];
        for (GUIPlotEquation plotEquation : plotEquations) {
            if (plotEquation == null) {
                continue;
            }
            double min = Math.max(xAxis.getLowerBound(), plotEquation.getMinX());
            double max = Math.min(xAxis.getUpperBound(), plotEquation.getMaxX());

            double delta = (max - min) / (nIncr + 1);
            double[] extras = plotEquation.getExtras();
            double[] ax = new double[1 + extras.length];
            for (int i = 0; i < nIncr; i++) {
                double xValue = min + (i + 1) * delta;
                double x = xAxis.getDisplayPosition(xValue);
                ax[0] = xValue;
                System.arraycopy(extras, 0, ax, 1, extras.length);
                double yValue = plotEquation.calculate(ax);
                yValue /= plotEquation.getScaleValue();
                double y = yAxis.getDisplayPosition(yValue);
                xValues[i] = x;
                yValues[i] = y;
            }
            gC.setStroke(plotEquation.getColor());
            gC.strokePolyline(xValues, yValues, nIncr);
        }
    }

    private ObservableList<DataSeries> loadChartData(String[] residues) throws IOException {
        ObservableList<DataSeries> data = FXCollections.observableArrayList();
        Path path = Paths.get(".", fileName);
        DataSeries series;
        simData.clear();
        double[] xValues = null;
        double[] yValues;
        try (BufferedReader fileReader = Files.newBufferedReader(path)) {
            String line = fileReader.readLine();
            double field = Double.parseDouble(line.trim());
            line = fileReader.readLine();
            double time = Double.parseDouble(line.trim());
            while (true) {
                line = fileReader.readLine();
                if (line == null) {
                    break;
                }
                String sline = line.trim();
                if (sline.isEmpty() || sline.startsWith("#")) {
                    continue;
                }
                String[] sfields = sline.split("\t");
                if (line.startsWith("vCPMG")) {
                    xValues = new double[sfields.length - 2];
                    for (int i = 2; i < sfields.length; i++) {
                        xValues[i - 2] = Double.parseDouble(sfields[i]);
                    }
                } else if (Arrays.stream(residues).anyMatch(e -> e.equals(sfields[0]))) {
                    series = new DataSeries();
                    data.add(series);

                    series.setName(sfields[0]);
                    double ref = Double.parseDouble(sfields[1]);
                    yValues = new double[sfields.length - 2];
                    for (int i = 2; i < sfields.length; i++) {
                        yValues[i - 2] = -Math.log(Double.parseDouble(sfields[i]) / ref) / time;
                        //-math.log(v/ref)/tCPMG
                    }
                    for (int i = 0; i < xValues.length; i++) {
                        double x = xValues[i];
                        double y = yValues[i];
                        XYValue dataPoint = new XYValue(x, y);
                        series.getData().add(dataPoint);
                    }
                }
            }
        }
        return data;
    }
    // Returns a list of maps that store data about each plot on cpmgfit

    public HashMap<String, Object>[] getPlottedData() {
        int numPlots = getNumPlots();
        int numFits = getNumEquations();

        if (numPlots == numFits) {
            // The fitted curve belongs to the raw data collected
            HashMap<String, Object>[] plottedData = new HashMap[numPlots];
            for (int i = 0; i < numPlots; i++) {
                HashMap<String, Object> plotDatum = new HashMap<>();
                plotDatum.put("rawData", returnLine(i));
                plotDatum.put("graphTitle", getSeriesName(i).replace(" ", "_")
                        .replace(".", "_"));
                plotDatum.put("fittedData", returnEquation(i));
                plottedData[i] = plotDatum;
            }
            return plottedData;

        }
        //Will break if numPlots != numFits
        return null;
    }

    // Returns a map of parameters about the graph
    public HashMap<String, Object> getGraphData() {
        HashMap<String, Object> graphData = new HashMap<>();

        String title = getTitle();
        String xLabel = xAxis.getLabel();
        String yLabel = yAxis.getLabel();
        double xMin = xAxis.lowerBoundProperty().getValue();
        double xMax = xAxis.upperBoundProperty().getValue();
        double yMin = yAxis.lowerBoundProperty().getValue();
        double yMax = yAxis.upperBoundProperty().getValue();

        List<Double> ranges = new ArrayList<>(Arrays.asList(xMin, xMax, yMin, yMax));
        graphData.put("title", title);
        graphData.put("xlabel", xLabel);
        graphData.put("ylabel", yLabel);
        graphData.put("ranges", ranges);
        graphData.put("colors", PlotData.colors);
        return graphData;
    }

    private double getScaledError(XYValue item) {
        Object extraValue = item.getExtraValue();
        double errorY = 0.0;
        if (extraValue instanceof ExperimentData.DataValue dataValue) {
            double error = dataValue.getError();
            double yValue = dataValue.getY();
            double errorY2 = yAxis.getDisplayPosition(yValue + error);
            double errorY1 = yAxis.getDisplayPosition(yValue - error);
            errorY = Math.abs(errorY2 - errorY1) / 2.0;
        } else if (extraValue instanceof Double error) {
            double yValue = (item.getYValue());
            double errorY2 = yAxis.getDisplayPosition(yValue + error);
            double errorY1 = yAxis.getDisplayPosition(yValue - error);
            errorY = Math.abs(errorY2 - errorY1) / 2.0;
        }
        return errorY;

    }

    public String getSeriesName(int iSeries) {
        ObservableList<DataSeries> data = getData();
        DataSeries series = data.get(iSeries);
        return series.getName();
    }

    public void dumpLine(int iSeries) {
        ArrayList<ArrayList<Double>> lineData = returnLine(iSeries);
        lineData.forEach((pointData -> {
            String outputLine = String.format("%14.4g %14.4g %14.4g", pointData.get(0), pointData.get(1), pointData.get(2));
            System.out.println(outputLine);
        }));
    }

    public ArrayList<ArrayList<Double>> returnLine(int iSeries) {
        ObservableList<DataSeries> data = getData();
        DataSeries series = data.get(iSeries);
        ArrayList<ArrayList<Double>> lineData = new ArrayList<>(series.getData().size());
        series.getData().forEach(value -> {
            Double x = value.getXValue();
            Double y = value.getYValue();
            Object extraValue = value.getExtraValue();
            if (extraValue instanceof ExperimentData.DataValue dataValue) {
                Double error = dataValue.getError();
                ArrayList<Double> pointData = new ArrayList<>(3);
                pointData.add(x);
                pointData.add(y);
                pointData.add(error);
                lineData.add(pointData);
            }
        });
        return lineData;
    }

    public List<List<Double>> returnEquation(int iLine) {
        PlotEquation plotEquation = plotEquations.get(iLine);
        double min = xAxis.getLowerBound();
        double max = xAxis.getUpperBound();
        int nIncr = 100;
        double delta = (max - min) / nIncr;
        List<List<Double>> equationData = new ArrayList<>(nIncr - 1);
        double[] extras = plotEquation.getExtras();
        double[] ax = new double[1 + extras.length];
        for (int i = 1; i < nIncr - 1; i++) {
            List<Double> pointData = new ArrayList<>(3);
            double xValue = min + i * delta;
            ax[0] = xValue;
            if (extras.length > 0) {
                ax[1] = extras[0];
            }
            if (extras.length - 1 >= 0) System.arraycopy(extras, 1, ax, 2, extras.length - 1);
            double yValue = plotEquation.calculate(ax);
            pointData.add(xValue);
            pointData.add(yValue);
            equationData.add(pointData);
        }
        return equationData;
    }

    public void saveEquationData() throws IOException {
        FileChooser fileChooser = new FileChooser();
        File file = fileChooser.showSaveDialog(null);
        if (file != null) {
            try (FileWriter fileWriter = new FileWriter(file)) {
                int numFits = getNumEquations();
                List<List<List<Double>>> allData = new ArrayList<>();
                for (int j = 0;j<numFits;j++) {
                    List<List<Double>> data = returnEquation(0);
                    allData.add(data);
                }
                for (int i = 0; i < allData.get(0).size(); i++) {
                    double x = allData.get(0).get(i).get(0);
                    StringBuilder outStr = new StringBuilder(String.format("%7f", x));
                    for (int j = 0;j<numFits;j++) {
                        double y = allData.get(j).get(i).get(1);
                        outStr.append(String.format(" %7f", y));
                    }
                    fileWriter.write(outStr + "\n");

                }
            }
        }

    }
}
