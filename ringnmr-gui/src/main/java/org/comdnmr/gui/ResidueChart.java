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

import java.util.*;

import javafx.collections.ListChangeListener;
import javafx.geometry.Orientation;
import javafx.scene.canvas.Canvas;
import javafx.scene.input.MouseEvent;
import javafx.scene.paint.Color;
import org.comdnmr.data.ExperimentSet;
import org.nmrfx.chemistry.relax.RelaxationSet;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.XYCanvasBarChart;
import org.nmrfx.chart.XYValue;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.ValueSet;
import org.nmrfx.graphicsio.GraphicsContextInterface;

/**
 *
 * @author brucejohnson
 */
public class ResidueChart extends XYCanvasBarChart {

    private static final Set<SelectionValue> selections = new LinkedHashSet<>();
    static List<String> mapNames = new ArrayList<>();
    public String currentSeriesName = "";
    ValueSet valueSet = null;
    Set<ResonanceSource> dynSources = new HashSet<>();

    public record SelectionValue(String seriesName, ResonanceSource resonanceSource) {}

    public static ResidueChart buildChart(Canvas canvas) {
        Axis xAxis = new Axis(Orientation.HORIZONTAL, 0, 100, 400, 100.0);
        Axis yAxis = new Axis(Orientation.VERTICAL, 0, 100, 100, 400);
        yAxis.setZeroIncluded(true);
        return new ResidueChart(canvas, xAxis, yAxis);
    }

    public ResidueChart(Canvas canvas, final Axis... AXIS) {
        super(canvas, AXIS);
        xAxis = AXIS[0];
        yAxis = AXIS[1];
        init();
    }

    public static List<ResonanceSource> getSelectedSources() {
        return selections.stream().map(sR ->sR.resonanceSource).distinct().toList();
    }

    public static Set<SelectionValue> getSelections() {
        return selections;
    }

    public static void clearSelections() {
        selections.clear();
    }

    public static void addSelections(List<SelectionValue> values) {
        selections.addAll(values);
    }

    private void init() {
        getData().addListener((ListChangeListener) (e -> seriesChanged()));
    }

    public String getExpSeries(String seriesName) {
        MoleculeBase mol = MoleculeFactory.getActive();

        if (seriesName.startsWith(mol.getName()) && seriesName.contains("_RING")) {
            String testSeries = seriesName.substring(mol.getName().length()+1, seriesName.indexOf("_RING"));
            if (ChartUtil.getResiduePropertyNames().contains(testSeries)) {
                ValueSet valueSet = ChartUtil.getResidueProperty(testSeries);
                seriesName = valueSet.name() + '|' +
                        "EXPAB" + "|" + "0:0:0" + "|" + "R" + "|" + "0";
            }
        }
        return seriesName;

    }
    public boolean mouseClicked(MouseEvent e) {
        Optional<Hit> hitOpt = pickChart(e.getX(), e.getY(), 5);
        boolean hitChart = false;
        if (hitOpt.isPresent()) {
            hitChart = true;
            Hit hit = hitOpt.get();
            boolean appendMode = e.isShiftDown();
            XYValue value = hit.getValue();
            var series = hit.getSeries();
            Object extraValue = value.getExtraValue();
            ResonanceSource resSource = null;
            if (extraValue instanceof ResonanceSource resonanceSource) {
                resSource = resonanceSource;
            }
            String seriesName = getExpSeries(series.getName());
            String statusMessage = seriesName + " " + resSource + " " + String.format("%.2f", value.getYValue());
            PyController.mainController.statusBar.setText(statusMessage);
            showInfo(seriesName, hit.getIndex(), resSource, appendMode);
        } else {
            Optional<ResonanceSource> intOpt = pickPresenceIndicators(e.getX(), e.getY());
            if (intOpt.isPresent()) {
                hitChart = true;
                boolean appendMode = e.isShiftDown();
                double f = e.getX() / xAxis.getWidth();
                int seriesIndex = (int) Math.floor(getData().size() * f);
                if ((seriesIndex >= 0) && (seriesIndex < getData().size())) {
                    String seriesName = getData().get(seriesIndex).getName();
                    showInfo(seriesName, seriesIndex, intOpt.get(), appendMode);
                }
            } else {
                clearSelections();
            }
        }
        return hitChart;
    }

    @Override
    public void annotate(GraphicsContextInterface gC) {
        if (dynSources != null && !dynSources.isEmpty()) {
            drawPresenceIndicators(gC);
        }
    }

    void seriesChanged() {
        clearSelections();
        drawChart();

    }

    public void setResProps(ValueSet valueSet) {
        this.valueSet = valueSet;
        dynSources = valueSet.resonanceSources();
    }

    void showInfo(String seriesName, int seriesIndex, ResonanceSource resSource, boolean appendMode) {
        if (!appendMode) {
            clearSelections();
            if (resSource != null) {
                selections.add(new SelectionValue(seriesName, resSource));
            }
            mapNames.clear();
        } else if (!getSelectedSources().contains(resSource)) {
            if (resSource != null) {
                selections.add(new SelectionValue(seriesName, resSource));
            }
        }
        currentSeriesName = seriesName;
        String[] seriesNameParts = seriesName.split("\\|");
        String mapName = seriesNameParts[0];
        mapNames.add(mapName);
        valueSet = ChartUtil.getResidueProperty(mapName);
        dynSources = valueSet.resonanceSources();
        showInfo(seriesName);
        drawChart();
    }

    void showInfo() {
        String[] seriesNameParts = currentSeriesName.split("\\|");
        String mapName = seriesNameParts[0];
        mapNames.clear();
        mapNames.add(mapName);
        valueSet = ChartUtil.getResidueProperty(mapName);
        dynSources = valueSet.resonanceSources();
        showInfo(currentSeriesName);
    }

    void showInfo(String seriesName) {
        PyController controller = PyController.mainController;
        PlotData xyCanvasChart = controller.xychart;
        String[] seriesNameParts = seriesName.split("\\|");
        if (seriesNameParts.length < 3) {
            return;
        }
        String equationName = seriesNameParts[1];
        String state = seriesNameParts[2];
        state = "*:" + state.substring(2);
//        System.out.println("series " + seriesName + " map " + mapName + " eqn " + equationName + " state " + state);
        SelectionValue[] resSources = new SelectionValue[selections.size()];
        selections.toArray(resSources);
        controller.chartInfo.valueSet = valueSet;
        controller.chartInfo.setResidues(resSources);
        controller.chartInfo.state = state;
        controller.chartInfo.mapName.clear();
        controller.chartInfo.mapName.addAll(mapNames);
        controller.chartInfo.equationName = equationName;

        if (valueSet instanceof ExperimentSet experimentSet) {
            if (!selections.isEmpty()) {
                controller.chartInfo.experimentalResult = experimentSet.getExperimentResult(resSources[0].resonanceSource);
            }
        } else if (valueSet instanceof RelaxationSet) {
        }
        controller.showInfo(controller.chartInfo, xyCanvasChart);
    }

    Optional<ResonanceSource> pickPresenceIndicators(double mouseX, double mouseY) {
        Optional<ResonanceSource> result = Optional.empty();
        if (dynSources != null) {
            for (var dynSource : dynSources) {
                int resNum = dynSource.getAtom().getResidueNumber();
                double x1 = xAxis.getDisplayPosition(resNum - 0.5) + 1;
                double x2 = xAxis.getDisplayPosition(resNum + 0.5) - 1;
                double y1 = yAxis.getYOrigin() - yAxis.getHeight() + 2;
                double y2 = yAxis.getYOrigin();
                if ((mouseX > x1) && (mouseX < x2)) {
                    if ((mouseY > y1) && (mouseY < y2)) {
                        result = Optional.of(dynSource);
                        break;
                    }
                }
            }
        }
        return result;
    }

    void drawPresenceIndicators(GraphicsContextInterface gC) {
        if (dynSources != null) {
            for (var dynSource : dynSources) {
                int resNum = dynSource.getAtom().getResidueNumber();
                double x1 = xAxis.getDisplayPosition(resNum - 0.5) + 1;
                double x2 = xAxis.getDisplayPosition(resNum + 0.5) - 1;
                double y1 = yAxis.getYOrigin() - yAxis.getHeight() + 2;
                double width = x2 - x1;
                double height = yAxis.getHeight() - 4;
                boolean hasDyn = selections.stream().map(sV -> sV.resonanceSource).anyMatch(sVres -> sVres == dynSource);
                if (hasDyn) {
                    gC.setFill(Color.ORANGE);
                } else {
                    gC.setFill(Color.LIGHTGRAY);
                }
                gC.fillRect(x1, y1, width, height);
            }
        }
    }
}
