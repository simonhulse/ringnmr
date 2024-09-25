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

import de.jensd.fx.glyphs.GlyphsDude;
import de.jensd.fx.glyphs.fontawesome.FontAwesomeIcon;
import javafx.animation.PauseTransition;
import javafx.application.Platform;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.geometry.Bounds;
import javafx.print.PrinterJob;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.SnapshotParameters;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.image.WritableImage;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.stage.StageStyle;
import javafx.util.Duration;
import org.comdnmr.data.*;
import org.comdnmr.eqnfit.*;
import org.comdnmr.fit.ResidueFitter;
import org.comdnmr.modelfree.CorrelationTime;
import org.comdnmr.modelfree.FitDeuteriumModel;
import org.comdnmr.modelfree.FitModel;
import org.comdnmr.modelfree.FitR1R2NOEModel;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.OrderParameterTool;
import org.comdnmr.util.CoMDOptions;
import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.util.ProcessingStatus;
import org.comdnmr.utils.NMRFxClient;
import org.controlsfx.control.StatusBar;
import org.controlsfx.dialog.ExceptionDialog;
import org.nmrfx.chart.Axis;
import org.nmrfx.chart.DataSeries;
import org.nmrfx.chart.XYEValue;
import org.nmrfx.chart.XYValue;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.InvalidMoleculeException;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.io.MoleculeIOException;
import org.nmrfx.chemistry.relax.*;
import org.nmrfx.chemistry.relax.RelaxTypes;
import org.nmrfx.console.ConsoleController;
import org.nmrfx.graphicsio.GraphicsIOException;
import org.nmrfx.graphicsio.SVGGraphicsContext;
import org.nmrfx.peaks.InvalidPeakException;
import org.nmrfx.star.ParseException;
import org.nmrfx.utils.GUIUtils;

import javax.imageio.ImageIO;
import javax.script.ScriptException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static org.comdnmr.gui.MainApp.preferencesController;
import static org.comdnmr.gui.MainApp.primaryStage;

public class PyController implements Initializable {

    public static PyController mainController;
    ResidueChart activeChart;
    Stage stage;
    List<ResidueChart> barCharts = new ArrayList<>();

    SSPainter ssPainter = null;

    @FXML
    StackPane stackPane;
    @FXML
    XYPlotDataPane chartPane;
    @FXML
    ScrollBar barScaleScrollBar;
    @FXML
    Slider barScaleSlider;

    @FXML
    Button nmrFxPeakButton;
    @FXML
    TableView<ExperimentData.DataValue> resInfoTable;
    @FXML
    TableView<ParValueInterface> parameterTable;
    @FXML
    ChoiceBox<String> equationChoice;
    @FXML
    TabPane parTabPane;
    @FXML
    Label aicLabel;
    @FXML
    Label aiccLabel;
    @FXML
    Label rmsLabel;
    @FXML
    Label rChiSqLabel;
    @FXML
    Label nLabel;

    @FXML
    ScrollPane chartBox;
    @FXML
    StatusBar statusBar;
    Circle statusCircle;
    @FXML
    Menu chartMenu;
    @FXML
    Menu experimentalDataAxisMenu;
    @FXML
    Menu moleculeDataAxisMenu;
    @FXML
    Menu groupsMenu;
    @FXML
    Menu orderParSetAxisMenu;
    @FXML
    BorderPane simPane;

    @FXML
    ChoiceBox<String> simChoice;
    @FXML
    CheckBox sliderGuessCheckBox;
    @FXML
    CheckBox calcErrorsCheckBox;

    @FXML
    TextField xLowerBoundTextField;
    @FXML
    TextField xUpperBoundTextField;
    @FXML
    TextField xTickTextField;
    @FXML
    TextField yLowerBoundTextField;
    @FXML
    TextField yUpperBoundTextField;
    @FXML
    TextField yTickTextField;
    @FXML
    CheckBox scalePlot;

    @FXML
    Button setBoundsButton;
    @FXML
    CheckBox autoscaleXCheckBox;
    @FXML
    CheckBox autoscaleYCheckBox;
    @FXML
    CheckBox zeroXCheckBox;
    @FXML
    CheckBox zeroYCheckBox;

    EquationControls simControls;

    @FXML
    TextField genDataNPtsTextField;
    @FXML
    TextField genDataXLBTextField;
    @FXML
    TextField genDataXUBTextField;
    @FXML
    TextField genDataSDevTextField;
    @FXML
    TextField genDataXValTextField;

    @FXML
    SplitPane splitPane;
    @FXML
    ChoiceBox<String> t1Choice;
    @FXML
    ChoiceBox<String> t2Choice;
    @FXML
    MenuButton relaxSetMenu;
    @FXML
    TextField r1MedianField;
    @FXML
    TextField r2MedianField;
    @FXML
    TextField tauCalcField;
    @FXML
    Slider sLambdaSlider;
    @FXML
    Label sLambdaLabel;
    @FXML
    Slider tauLambdaSlider;
    @FXML
    Label tauLambdaLabel;
    @FXML
    CheckBox fitJCheckBox;
    @FXML
    ChoiceBox<FitModel.BootstrapMode> bootStrapChoice;
    @FXML
    CheckBox lambdaCheckBox;
    @FXML
    Slider tauFractionSlider;
    @FXML
    Label tauFractionLabel;
    @FXML
    Slider t2LimitSlider;
    @FXML
    Label t2LimitLabel;
    @FXML
    Slider nReplicatesSlider;
    @FXML
    Label bootstrapNLabel;
    @FXML
    HBox modelBox;
    List<CheckBox> modelCheckBoxes = new ArrayList<>();

    BootstrapSamplePlots bootstrapSamplePlots = null;
    InputDataInterface inputDataInterface = null;
    NMRFxClient cl;

    ResidueFitter residueFitter;
    FitModel modelFitter;
    List<ResonanceSource> fittingResidues = new ArrayList<>();
    boolean simulate = true;
    ChartInfo chartInfo = new ChartInfo();

    FitResult fitResult;
    PlotData xychart;
    Canvas barPlotCanvas = new Canvas();
    Optional<Double> barChartXMin = Optional.empty();
    Optional<Double> barChartXMax = Optional.empty();

    @FXML
    ToolBar navigatorToolBar = new ToolBar();

    static Random rand = new Random();
    File initialDir = null;
    SeriesComparator seriesComparator = new SeriesComparator();

    Function<String, String> nmrfxFunction;

    @FXML
    private void pyAction(ActionEvent event) {
        Node node = (Node) event.getSource();
        String id = node.getId();
        try {
            MainApp.engine.eval("print('howdy " + id + "')");
            MainApp.engine.eval("onAction('" + id + "')");
        } catch (ScriptException ignored) {

        }
        //MainApp.interpreter.exec("onAction(" + node + ")");
    }

    public void setNMRFxFunction(Function<String, String> nmrfxFunction) {
        this.nmrfxFunction = nmrfxFunction;
        nmrFxPeakButton.setDisable(false);
    }

    @FXML
    public void displayEquation() {
        try {
            simControls.simSliderAction("");
        } catch (NullPointerException npE) {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Residue must be selected in display equation");
            alert.showAndWait();
        }

    }

    public static PyController create(Stage stage) {
        FXMLLoader loader = new FXMLLoader(PyController.class.getResource("/fxml/RINGScene.fxml"));
        PyController controller = null;

        if (stage == null) {
            stage = new Stage(StageStyle.DECORATED);
        }

        try {
            Scene scene = new Scene(loader.load());
            stage.setScene(scene);
            scene.getStylesheets().add("/styles/Styles.css");

            controller = loader.getController();
            controller.stage = stage;
            //controllers.add(controller);
            stage.setTitle("RING NMR Dynamics");
            stage.show();
        } catch (IOException ioE) {

            ioE.printStackTrace();
            System.out.println(ioE.getMessage());
        }
        return controller;
    }

    @Override
    public void initialize(URL url, ResourceBundle rb) {
        mainController = this;
        makeAxisMenu();
        groupsMenu.showingProperty().addListener(e -> addOrderParSetsToAxisMenu());

        if (getFittingMode().equals("cpmg")) {
            simControls = new CPMGControls();
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1100.0");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("65.0");
            xTickTextField.setText("100.0");
            yTickTextField.setText("5.0");
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(xval);
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getFittingMode().equals("cest")) {
            simControls = new CESTControls();
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (hasExperimentSet()) {
                if (getCurrentExperimentSet().getExperimentData() != null) {
                    getCurrentExperimentSet().getExperimentData().
                            stream().findFirst().ifPresent(expData -> {
                                double[] xVals = ((DoubleArrayExperiment) expData).getXVals();
                                xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                                xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
                            });
                }
            }

            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.0");
            xTickTextField.setText("1.0");
            yTickTextField.setText("0.25");
            genDataNPtsTextField.setDisable(false);
            genDataNPtsTextField.setText("100");
            genDataXLBTextField.setDisable(false);
            genDataXLBTextField.setText("-8.0");
            genDataXUBTextField.setDisable(false);
            genDataXUBTextField.setText("8.0");
            genDataXValTextField.setDisable(true);
            genDataXValTextField.setText("");
        } else if (getFittingMode().equals("exp")) {
            simControls = new ExpControls();
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("100.0");
            xTickTextField.setText("0.25");
            yTickTextField.setText("10.0");
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(xval);
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getFittingMode().equals("noe")) {
            simControls = new NOEControls();
            xLowerBoundTextField.setText("-0.5");
            xUpperBoundTextField.setText("1.5");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.25");
            xTickTextField.setText("0.25");
            yTickTextField.setText("0.25");
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(xval);
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getFittingMode().equals("r1rho")) {
            simControls = new R1RhoControls();
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (hasExperimentSet()) {
                if (getCurrentExperimentSet().getExperimentData() != null) {
                    getCurrentExperimentSet().getExperimentData().
                            stream().findFirst().ifPresent(expData -> {
                                double[] xVals = ((DoubleArrayExperiment) expData).getXVals();
                                xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                                xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
                            });
                }
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("50.0");
            xTickTextField.setText("1.0");
            yTickTextField.setText("0.25");
            genDataNPtsTextField.setDisable(false);
            genDataNPtsTextField.setText("100");
            genDataXLBTextField.setDisable(false);
            genDataXLBTextField.setText("-8.0");
            genDataXUBTextField.setDisable(false);
            genDataXUBTextField.setText("8.0");
            genDataXValTextField.setDisable(true);
            genDataXValTextField.setText("");
        } else {
            System.out.println("Error: no fitting mode selected.");
        }
        VBox vBox = simControls.makeControls(mainController);
        simPane.centerProperty().set(vBox);
        CoMDOptions options = new CoMDOptions(true);
        residueFitter = new ResidueFitter(options, this::updateFitProgress, this::updateStatus);
        statusCircle = new Circle(10);
        statusBar.getLeftItems().add(statusCircle);
        updateEquationChoices(getFittingMode());
        equationChoice.valueProperty().addListener(e -> equationAction());

        simChoice.getItems().add("Simulate CPMG");
        simChoice.getItems().add("Simulate EXP");
        simChoice.getItems().add("Simulate CEST");
        simChoice.getItems().add("Simulate R1Rho");
        simChoice.setValue("Simulate CPMG");
        simChoice.valueProperty().addListener(s -> simAction());

        splitPane.setDividerPositions(0.4, 0.7);
        relaxSetMenu.showingProperty().addListener(e -> makeRelaxSetMenu());

        setBoundsButton.setOnAction(event -> setBounds());

        initResidueNavigator();
        calcErrorsCheckBox.selectedProperty().addListener(e -> FitFunction.setCalcError(calcErrorsCheckBox.isSelected()));
        calcErrorsCheckBox.setSelected(true);
        autoscaleXCheckBox.setSelected(true);
        autoscaleXCheckBox.selectedProperty().addListener(e -> autoscaleX(autoscaleXCheckBox.isSelected()));
        autoscaleYCheckBox.setSelected(true);
        autoscaleYCheckBox.selectedProperty().addListener(e -> autoscaleY(autoscaleYCheckBox.isSelected()));
        zeroXCheckBox.selectedProperty().addListener(e -> includeXZero(zeroXCheckBox.isSelected()));
        zeroXCheckBox.setSelected(false);
        zeroYCheckBox.selectedProperty().addListener(e -> includeYZero(zeroYCheckBox.isSelected()));
        zeroYCheckBox.setSelected(false);
        nmrFxPeakButton.setDisable(true);
        nmrFxPeakButton.setOnAction(this::nmrFxMessage);
        xychart = (PlotData) chartPane.getChart();
        xychart.getXAxis().setAutoRanging(true);
        xychart.getYAxis().setAutoRanging(true);
        chartPane.widthProperty().addListener(e -> resizeXYPlotCanvas());
        chartPane.heightProperty().addListener(e -> resizeXYPlotCanvas());
        chartBox.widthProperty().addListener(e -> resizeBarPlotCanvas());
        chartBox.heightProperty().addListener(e -> resizeBarPlotCanvas());
        chartBox.setContent(barPlotCanvas);
        addChart();
        barPlotCanvas.setOnMouseClicked(this::mouseClickedOnBarCanvas);
//        mainController.setOnHidden(e -> Platform.exit());
        PauseTransition logoTransition = new PauseTransition(Duration.seconds(5));
        logoTransition.setOnFinished(e -> removeLogo());
        logoTransition.play();
        chartMenu.setOnShowing(e -> makeAxisMenu());
        barScaleSlider.valueProperty().addListener(e -> resizeBarPlotCanvas());
        barScaleScrollBar.valueProperty().addListener(e -> resizeBarPlotCanvas());

        tauFractionSlider.setMin(0.0);
        tauFractionSlider.setMax(0.9);
        tauFractionSlider.setValue(0.1);
        tauFractionLabel.setText("0.1");
        tauFractionSlider.setBlockIncrement(0.1);
        tauFractionSlider.setMinorTickCount(4);
        tauFractionSlider.setMajorTickUnit(0.1);
        tauFractionSlider.setShowTickMarks(true);
        tauFractionSlider.setShowTickLabels(true);
        tauFractionSlider.valueProperty().addListener(v
                -> tauFractionLabel.setText(String.format("%.2f", tauFractionSlider.getValue())));


        sLambdaSlider.setMin(0.0);
        sLambdaSlider.setMax(2.0);
        sLambdaSlider.setValue(0.0);
        sLambdaLabel.setText("0.0");
        sLambdaSlider.setBlockIncrement(0.05);
        sLambdaSlider.setMajorTickUnit(0.5);
        sLambdaSlider.setMinorTickCount(4);
        sLambdaSlider.setShowTickMarks(true);
        sLambdaSlider.setShowTickLabels(true);
        sLambdaSlider.valueProperty().addListener(v
                -> sLambdaLabel.setText(String.format("%.1f", sLambdaSlider.getValue())));

        tauLambdaSlider.setMin(0.0);
        tauLambdaSlider.setMax(2.0);
        tauLambdaSlider.setValue(0.0);
        tauLambdaLabel.setText("0.0");
        tauLambdaSlider.setBlockIncrement(0.05);
        tauLambdaSlider.setMajorTickUnit(0.5);
        tauLambdaSlider.setMinorTickCount(4);
        tauLambdaSlider.setShowTickMarks(true);
        tauLambdaSlider.setShowTickLabels(true);
        tauLambdaSlider.valueProperty().addListener(v
                -> tauLambdaLabel.setText(String.format("%.1f", tauLambdaSlider.getValue())));


        t2LimitSlider.setMin(0.0);
        t2LimitSlider.setMax(100.0);
        t2LimitSlider.setValue(0.0);
        t2LimitLabel.setText("0.0");
        t2LimitSlider.setBlockIncrement(1.0);
        t2LimitSlider.setMajorTickUnit(25.0);
        t2LimitSlider.setMinorTickCount(4);
        t2LimitSlider.setShowTickMarks(true);
        t2LimitSlider.setShowTickLabels(true);
        t2LimitSlider.valueProperty().addListener(v
                -> t2LimitLabel.setText(String.format("%.1f", t2LimitSlider.getValue())));

        nReplicatesSlider.setMin(0.0);
        nReplicatesSlider.setMax(1000);
        nReplicatesSlider.setValue(0.0);
        bootstrapNLabel.setText("0");
        nReplicatesSlider.setBlockIncrement(25);
        nReplicatesSlider.setMajorTickUnit(50);
        nReplicatesSlider.setMinorTickCount(4);
        nReplicatesSlider.setShowTickMarks(true);
        nReplicatesSlider.setShowTickLabels(true);
        nReplicatesSlider.valueProperty().addListener(v
                -> bootstrapNLabel.setText(String.format("%d", (int) nReplicatesSlider.getValue())));
        bootStrapChoice.getItems().addAll(FitModel.BootstrapMode.values());
        bootStrapChoice.setValue(FitModel.BootstrapMode.PARAMETRIC);
        bootStrapChoice.valueProperty().addListener(e -> bootStrapChanged(bootStrapChoice.getValue()));
        lambdaCheckBox.selectedProperty().addListener(e -> lambdaChanged(lambdaCheckBox.isSelected()));

        modelBox.setSpacing(10);
        String[] modelNames = {"1", "1f", "1s", "2f", "2s", "1sf", "2sf"};
        for (var modelName : modelNames) {
            var checkBox = new CheckBox(modelName);
            checkBox.setMinWidth(40.0);
            modelBox.getChildren().add(checkBox);
            modelCheckBoxes.add(checkBox);
        }
    }

    public Stage getStage() {
        return stage;
    }

    void removeLogo() {
        stackPane.getChildren().remove(1);
    }

    void resizeXYPlotCanvas() {
        Canvas canvas = xychart.getCanvas();
        GraphicsContext gC = canvas.getGraphicsContext2D();
        canvas.setWidth(chartPane.getWidth());
        canvas.setHeight(chartPane.getHeight());
        gC.clearRect(0, 0, canvas.getWidth(), canvas.getHeight());
        xychart.setWidth(canvas.getWidth());
        xychart.setHeight(canvas.getHeight());
    }

    void filterSeries() {
        for (ResidueChart residueChart : barCharts) {
            for (DataSeries series : residueChart.getData()) {
                var copyOfValues = List.copyOf(series.getValues());
                series.clear();
                for (var v : copyOfValues) {
                    Object obj = v.getExtraValue();
                    if (obj instanceof ResonanceSource resonanceSource) {
                        if (!resonanceSource.deleted()) {
                            series.add(v);
                        }
                    }
                }
            }
        }
    }

    void resizeBarPlotCanvas() {
        double width = chartBox.getWidth();
        double height = chartBox.getHeight();
        double barScale = barScaleSlider.getValue();
        double barValue = barScaleScrollBar.getValue();
        barScaleScrollBar.setVisible(true);
        double scrollMin = barScaleScrollBar.getMin();
        double scrollMax = barScaleScrollBar.getMax();
        double visAmount = (scrollMax - scrollMin) / barScale;
        barScaleScrollBar.setVisibleAmount(visAmount);
        barPlotCanvas.setWidth(width);
        barPlotCanvas.setHeight(height);
        double visFraction = 1.0 / barScale;
        double barCenter = barValue / barScaleScrollBar.getMax() * (1.0 - visFraction);
        double barDiv = 1.0 / barScale;
        double fMax = barCenter + barDiv;
        GraphicsContext gC = barPlotCanvas.getGraphicsContext2D();
        gC.clearRect(0, 0, barPlotCanvas.getWidth(), barPlotCanvas.getHeight());
        if (ssPainter != null) {
            height -= ssPainter.getHeight();
        }
        double xMin = Double.MAX_VALUE;
        double xMax = Double.NEGATIVE_INFINITY;
        xMin = Math.floor((ChartUtil.minRes - 2) / 5.0) * 5.0;
        xMax = Math.ceil((ChartUtil.maxRes + 2) / 5.0) * 5.0;
        barChartXMin = Optional.of(xMin);
        barChartXMax = Optional.of(xMax);

        System.out.println("xx " + xMin + " " + xMax);
        double limitMin = barCenter * (xMax - xMin) + xMin;
        double limitMax = fMax * (xMax - xMin) + xMin;
        for (ResidueChart residueChart : barCharts) {
            for (DataSeries series : residueChart.getData()) {
                series.setLimits(limitMin, limitMax);
            }
        }
        int iChart = 0;
        double noAxisHeight = 0;
        double withAxisHeight = 0;
        for (ResidueChart residueChart : barCharts) {
            residueChart.xAxis.setLowerBound(limitMin);
            residueChart.xAxis.setUpperBound(limitMax);

            residueChart.setWidth(width);
            residueChart.autoScale(false);
            iChart++;
            boolean xVisible = iChart == barCharts.size();
            residueChart.getXAxis().setLabelVisible(xVisible);
            residueChart.getXAxis().setTickLabelsVisible(xVisible);
            if (xVisible) {
                withAxisHeight = residueChart.getXAxis().getBorderSize();
            } else {
                noAxisHeight = residueChart.getXAxis().getBorderSize();
            }
        }
        double extraHeight = withAxisHeight - noAxisHeight;
        double chartHeight = (height - extraHeight) / barCharts.size();
        double yPos = 0.0;
        for (ResidueChart residueChart : barCharts) {
            if (residueChart == barCharts.get(barCharts.size() - 1)) {
                residueChart.setHeight(chartHeight + extraHeight);
            } else {
                residueChart.setHeight(chartHeight);
            }
            residueChart.setYPos(yPos);
            residueChart.drawChart();
            yPos += chartHeight;
        }
        refreshResidueCharts();
    }

    void setResidueChartBorder() {
        double maxLeftBorder = 0.0;
        for (ResidueChart residueChart : barCharts) {
            double[] borders = residueChart.getMinBorders();
            maxLeftBorder = Math.max(borders[0], maxLeftBorder);
        }
        for (ResidueChart residueChart : barCharts) {
            residueChart.setMinLeftBorder(maxLeftBorder);
        }
    }

    void refreshResidueCharts() {
        setResidueChartBorder();
        GraphicsContext gC = barPlotCanvas.getGraphicsContext2D();
        gC.clearRect(0, 0, barPlotCanvas.getWidth(), barPlotCanvas.getHeight());
        for (ResidueChart residueChart : barCharts) {
            residueChart.drawChart();
        }
        if (ssPainter != null) {
            double ssHeight = ssPainter.getHeight();
            Axis axis = activeChart.xAxis;
            ssPainter.paintSS(axis.getXOrigin(), barPlotCanvas.getHeight() - ssHeight,
                    axis.getWidth(),
                    axis.getLowerBound(), axis.getUpperBound());
        }
    }

    void mouseClickedOnBarCanvas(MouseEvent e) {
        double x = e.getX();
        double y = e.getY();
        for (ResidueChart residueChart : barCharts) {
            Axis xAxis = residueChart.xAxis;
            Axis yAxis = residueChart.yAxis;
            if (residueChart.mouseClicked(e)) {
                activeChart = residueChart;
                System.out.println("sources " + residueChart.getSelectedSources());
                System.out.println("chartinfo " + chartInfo);

                break;
            } else if ((x > xAxis.getXOrigin()) && (x < xAxis.getXOrigin() + xAxis.getWidth())) {
                if ((y < yAxis.getYOrigin()) && (y > xAxis.getYOrigin() - yAxis.getHeight())) {
                    activeChart = residueChart;
                    break;
                }
            }
        }
        refreshResidueCharts();
    }

    public void setControls() {
        boolean update = false;
        if (getFittingMode().equals("cpmg") && !(simControls instanceof CPMGControls)) {
            simControls = new CPMGControls();
            update = true;
        } else if (getFittingMode().equals("exp") && !(simControls instanceof ExpControls)) {
            simControls = new ExpControls();
            update = true;
        } else if (getFittingMode().equals("noe") && !(simControls instanceof NOEControls)) {
            simControls = new NOEControls();
            update = true;
        } else if (getFittingMode().equals("cest") && !(simControls instanceof CESTControls)) {
            simControls = new CESTControls();
            update = true;
        } else if (getFittingMode().equals("r1rho") && !(simControls instanceof R1RhoControls)) {
            simControls = new R1RhoControls();
            update = true;
        }
        if (update) {
            VBox vBox = simControls.makeControls(mainController);
            simPane.centerProperty().set(vBox);
            updateEquationChoices(getFittingMode());
        }
        if (hasExperimentSet()) {
            getCurrentExperimentSet().getExperimentData().stream().findFirst().ifPresent(e -> simControls.setNucleus(e.getNucleusName()));
            updateXYChartLabels();
        }
    }

    public void setSimControls() {
        boolean update = false;
        if (getSimMode().equals("cpmg") && !(simControls instanceof CPMGControls)) {
            simControls = new CPMGControls();
            update = true;
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(xval);
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getSimMode().equals("exp") && !(simControls instanceof ExpControls)) {
            simControls = new ExpControls();
            update = true;
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(xval);
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getSimMode().equals("noe") && !(simControls instanceof NOEControls)) {
            simControls = new NOEControls();
            update = true;
            genDataNPtsTextField.setDisable(true);
            genDataXLBTextField.setDisable(true);
            genDataXUBTextField.setDisable(true);
            genDataXValTextField.setDisable(false);
            StringBuilder sBuilder = new StringBuilder();
            for (double xval : getFitter().getSimXDefaults()) {
                sBuilder.append(xval);
                sBuilder.append(" ");
            }
            genDataXValTextField.setText(sBuilder.toString());
        } else if (getSimMode().equals("cest") && !(simControls instanceof CESTControls)) {
            simControls = new CESTControls();
            ((CESTControls) simControls).updateDeltaLimits();
            update = true;
            genDataNPtsTextField.setDisable(false);
            genDataNPtsTextField.setText("100");
            genDataXLBTextField.setDisable(false);
            genDataXLBTextField.setText("-8.0");
            genDataXUBTextField.setDisable(false);
            genDataXUBTextField.setText("8.0");
            genDataXValTextField.setDisable(true);
            genDataXValTextField.setText("");
        } else if (getSimMode().equals("r1rho") && !(simControls instanceof R1RhoControls)) {
            simControls = new R1RhoControls();
            ((R1RhoControls) simControls).updateDeltaLimits();
            update = true;
            genDataNPtsTextField.setDisable(false);
            genDataNPtsTextField.setText("100");
            genDataXLBTextField.setDisable(false);
            genDataXLBTextField.setText("-8.0");
            genDataXUBTextField.setDisable(false);
            genDataXUBTextField.setText("8.0");
            genDataXValTextField.setDisable(true);
            genDataXValTextField.setText("");
        }
        if (update) {
            updateEquationChoices(getSimMode());
            VBox vBox = simControls.makeControls(mainController);
            simPane.centerProperty().set(vBox);
        }
    }

    void updateEquationChoices() {
        updateEquationChoices(getFittingMode());
    }

    void updateEquationChoices(String mode) {
        switch (mode) {
            case "cpmg":
                simControls.updateEquations(equationChoice, CPMGFitter.getEquationNames());
                break;
            case "exp":
                simControls.updateEquations(equationChoice, ExpFitter.getEquationNames());
                break;
            case "noe":
                simControls.updateEquations(equationChoice, NOEFit.getEquationNames());
                break;
            case "cest":
                if (!CESTFitter.getEquationNames().isEmpty()) {
                    simControls.updateEquations(equationChoice, CESTFitter.getEquationNames());
                } else {
                    Alert alert = new Alert(Alert.AlertType.ERROR, "At least one CEST equation must be selected in Preferences.");
                    alert.showAndWait();
                }
                break;
            case "r1rho":
                if (!R1RhoFitter.getEquationNames().isEmpty()) {
                    simControls.updateEquations(equationChoice, R1RhoFitter.getEquationNames());
                } else {
                    Alert alert = new Alert(Alert.AlertType.ERROR, "At least one R1RHO equation must be selected in Preferences.");
                    alert.showAndWait();
                }
                break;
        }
    }

    void initResidueNavigator() {

        String iconSize = "12px";
        String fontSize = "7pt";
        ArrayList<Button> buttons = new ArrayList<>();
        Button bButton;

        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.FAST_BACKWARD, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(this::firstResidue);
        buttons.add(bButton);
        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.BACKWARD, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(this::previousResidue);
        buttons.add(bButton);
        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.FORWARD, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(this::nextResidue);
        buttons.add(bButton);
        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.FAST_FORWARD, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(this::lastResidue);
        buttons.add(bButton);
        bButton = GlyphsDude.createIconButton(FontAwesomeIcon.REMOVE, "", iconSize, fontSize, ContentDisplay.GRAPHIC_ONLY);
        bButton.setOnAction(this::removeItem);
        buttons.add(bButton);

        navigatorToolBar.getItems().addAll(buttons);
    }

    private void incrResidue(int delta) {
        var currentSet = getCurrentExperimentSet();
        if (currentSet == null) {
            return;
        }
        List<ExperimentResult> resInfo = currentSet.getExperimentResults();
        List<ResonanceSource> resSources = new ArrayList<>();
        for (ExperimentResult experimentResult : resInfo) {
            resSources.add(experimentResult.getResonanceSource());
        }
        Collections.sort(resSources);
        if (chartInfo.hasResidue()) {
            int resIndex = resSources.indexOf(chartInfo.getSource());
            resIndex += delta;
            if (resIndex <= 0) {
                resIndex = 0;
            } else if (resIndex >= resSources.size()) {
                resIndex = resSources.size() - 1;
            }
            var resSource = resSources.get(resIndex);
            ResidueChart chart = getActiveChart();
            chart.showInfo(chart.currentSeriesName, 0, resSource, false);
            String statusMessage = chart.currentSeriesName + " " + resSource;
            mainController.statusBar.setText(statusMessage);

        }
    }

    public void previousResidue(ActionEvent event) {
        incrResidue(-1);
    }

    public void firstResidue(ActionEvent event) {
        incrResidue(-10000);
    }

    public void nextResidue(ActionEvent event) {
        incrResidue(1);
    }

    public void lastResidue(ActionEvent event) {
        incrResidue(10000);
    }

    void removeItem(ActionEvent event) {
        if (chartInfo != null) {
            for (var resonanceSource : chartInfo.getResidues()) {
                resonanceSource.deleted(!resonanceSource.deleted());
            }
            filterSeries();
            refreshResidueCharts();
        }
    }

    public void loadParameterFile() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open YAML File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("Yaml File", "*.yaml", "*.yml"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        if (file != null) {
            if (activeChart != null) {
                clearChart();
                simulate = false;
                fitResult = null;
            }
            chartInfo.clear();
            ChartUtil.loadParameters(file.toString());
        }
        clearSecondaryStructure();
    }

    @FXML
    public void loadPDBFile() throws MoleculeIOException, ParseException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open PDB File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("PDB File", "*.pdb"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        String type = "pdb";
        if (file != null) {
            ChartUtil.loadMoleculeFile(file.toString(), type);
        }
    }

    @FXML
    public void loadSTARFile() throws MoleculeIOException, ParseException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open STAR File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("STAR File", "*.str"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        String type = "star";
        if (file != null) {
            ChartUtil.loadMoleculeFile(file.toString(), type);
        }
    }

    @FXML
    public void loadNEFFile() throws MoleculeIOException, ParseException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open NEF File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("NEF File", "*.nef"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        String type = "nef";
        if (file != null) {
            ChartUtil.loadMoleculeFile(file.toString(), type);
        }
    }

    @FXML
    public void loadCIFFile() throws MoleculeIOException, ParseException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open CIF File");
        fileChooser.getExtensionFilters().addAll(new FileChooser.ExtensionFilter("CIF File", "*.cif"));
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        String type = "cif";
        if (file != null) {
            ChartUtil.loadMoleculeFile(file.toString(), type);
        }
    }

    @FXML
    public void loadRelaxValues() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Open Relaxation Values File");
        Stage stage = MainApp.primaryStage;
        File file = fileChooser.showOpenDialog(stage);
        if (file != null) {
            try {
                DataIO.loadRelaxationTextFile(file);
            } catch (IOException ioE) {
                ExceptionDialog dialog = new ExceptionDialog(ioE);
                dialog.showAndWait();
                return;
            }
            addMoleculeDataToAxisMenu();
            showT1T2NoeData();
        }
    }

    @FXML
    public void loadSecondaryStructure() {
        List<SecondaryStructure> ssValues = SecondaryStructure.loadFromFile();
        if (!ssValues.isEmpty()) {
            ssPainter = new SSPainter(barPlotCanvas, ssValues);
        }
        resizeBarPlotCanvas();
    }

    @FXML
    public void clearSecondaryStructure() {
        if (ssPainter != null) {
            ssPainter.secondaryStructures.clear();
            ssPainter = null;
            resizeBarPlotCanvas();
        }
    }

    @FXML
    public void inputParameters() {
        if (inputDataInterface == null) {
            inputDataInterface = new InputDataInterface(this);
        }
        inputDataInterface.inputParameters();
    }

    @FXML
    public void loadFromPeakLists() {
        if (inputDataInterface == null) {
            inputDataInterface = new InputDataInterface(this);
        }
        inputDataInterface.createPeakListInterface();
    }

    @FXML
    public void startServer() {
        String tempDir = System.getProperty("java.io.tmpdir");
        String userName = System.getProperty("user.name");
        Path path = FileSystems.getDefault().getPath(tempDir, "NMRFx_" + userName + "_port.txt");
        File f = path.toFile(); //new File(tempDir + "/NMRFx_" + userName + "_port.txt");
        int port = 8021;
        if (!f.exists() && !f.isDirectory()) {
            Alert alert = new Alert(Alert.AlertType.WARNING, "NMRFx Server port file " + "\n" + f + "\ndoesn't exist. \nCreate NMRFx server first in NMRFxProcessor GUI.");
            alert.showAndWait();
        } else {
            try {
                BufferedReader reader = new BufferedReader(new FileReader(f));
                String text;
                while ((text = reader.readLine()) != null) {
                    port = Integer.parseInt(text);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            cl = NMRFxClient.makeClient(port);
            nmrFxPeakButton.setDisable(false);
        }
    }

    public NMRFxClient getClient() {
        return cl;
    }

    @FXML
    void nmrFxMessage(ActionEvent e) {
        String peakString = getPeakNumFromTable();
        System.out.println("show peak " + peakString);
        nmrfxFunction.apply("nw.showPeak " + peakString);
    }

    public void updateXYChartLabels() {
        if ((simControls instanceof CPMGControls)) {
            xychart.setNames("CPMG", "\u03BD (cpmg)", "R2 (\u03BD)", "0");
            xychart.setBounds(0.0, 1100.0, 0.0, 65.0, 100.0, 5.0);
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1100.0");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("65.0");
            xTickTextField.setText("100.0");
            yTickTextField.setText("5.0");
        } else if ((simControls instanceof ExpControls)) {
            xychart.setNames("Exp", "Time (s)", "Intensity", "0");
            xychart.setBounds(0.0, 1.25, 0.0, 100, 0.25, 10.0);
            xLowerBoundTextField.setText("0.0");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("100.0");
            xTickTextField.setText("0.25");
            yTickTextField.setText("10.0");
        } else if ((simControls instanceof NOEControls)) {
            xychart.setNames("NOE", "On/Off", "Intensity", "0");
            xychart.setBounds(-0.25, 1.25, 0.0, 1.25, 0.25, 10.0);
            xLowerBoundTextField.setText("00.25");
            xUpperBoundTextField.setText("1.25");
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.25");
            xTickTextField.setText("0.25");
            yTickTextField.setText("0.25");
        } else if ((simControls instanceof CESTControls)) {
            xychart.setNames("CEST", "Offset (PPM)", "I(t)/I(0)", "20");
            xychart.setBounds(-20, 20, 0.0, 1.0, 2.0, 0.25);
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (hasExperimentSet()) {
                getCurrentExperimentSet().getExperimentData().
                        stream().findFirst().ifPresent(e -> {
                            double[] xVals = ((DoubleArrayExperiment) e).getXVals();
                            if (xVals != null) {
                                xychart.setBounds(Math.floor(xVals[1] / 2) * 2, Math.ceil(xVals[xVals.length - 1] / 2) * 2, 0.0, 1.0, 1.0, 0.25);
                                xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                                xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
                            }
                        });
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("1.0");
            xTickTextField.setText("2.0");
            yTickTextField.setText("0.25");
        } else if ((simControls instanceof R1RhoControls)) {
            xychart.setNames("R1Rho", "Offset (PPM)", "R1Rho", "20");
            xychart.setBounds(-20, 20, 0.0, 50.0, 2.0, 5.0);
            xLowerBoundTextField.setText("-20.0");
            xUpperBoundTextField.setText("20.0");
            if (hasExperimentSet()) {
                getCurrentExperimentSet().getExperimentData().
                        stream().findFirst().ifPresent(e -> {
                            double[] xVals = ((DoubleArrayExperiment) e).getXVals();
                            if (xVals != null) {
                                xychart.setBounds(Math.floor(xVals[1] / 2) * 2, Math.ceil(xVals[xVals.length - 1] / 2) * 2, 0.0, 1.0, 1.0, 0.25);
                                xLowerBoundTextField.setText(String.valueOf(Math.floor(xVals[1] / 2) * 2));
                                xUpperBoundTextField.setText(String.valueOf(Math.ceil(xVals[xVals.length - 1] / 2) * 2));
                            }
                        });
            }
            yLowerBoundTextField.setText("0.0");
            yUpperBoundTextField.setText("50.0");
            xTickTextField.setText("2.0");
            yTickTextField.setText("5.0");
        }
    }

    public void setBounds() {
        try {
            double xLB = Double.parseDouble(xLowerBoundTextField.getText());
            double xUB = Double.parseDouble(xUpperBoundTextField.getText());
            double xTick = Double.parseDouble(xTickTextField.getText());
            double yLB = Double.parseDouble(yLowerBoundTextField.getText());
            double yUB = Double.parseDouble(yUpperBoundTextField.getText());
            double yTick = Double.parseDouble(yTickTextField.getText());
            if (xLB < xUB && yLB < yUB) {
                xychart.setBounds(xLB, xUB, yLB, yUB, xTick, yTick);
            } else {
                Alert alert = new Alert(Alert.AlertType.ERROR);
                alert.setContentText("Error: Upper Bound must be > Lower Bound.");
                alert.showAndWait();
            }

        } catch (NumberFormatException nfEsb) {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Bound and Interval values must be provided.");
            alert.showAndWait();
        }

    }

    // temporary addition till bug in xychart autoscale (only uses one series) is fixed
    public double[] calcAutoScale() {
        if (!xychart.getData().isEmpty()) {
            double xMax = Double.NEGATIVE_INFINITY;
            double xMin = Double.MAX_VALUE;
            double yMax = Double.NEGATIVE_INFINITY;
            double yMin = Double.MAX_VALUE;
            boolean ok = false;
            for (DataSeries dataSeries : xychart.getData()) {
                if (!dataSeries.isEmpty()) {
                    ok = true;
                    xMin = Math.min(xMin, dataSeries.getMinX());
                    xMax = Math.max(xMax, dataSeries.getMaxX());
                    yMin = Math.min(yMin, dataSeries.getMinY());
                    yMax = Math.max(yMax, dataSeries.getMaxY());
                }
            }

            if (ok) {
                double[] bounds = {xMin, xMax, yMin, yMax};
                return bounds;
            }
        }
        return null;
    }

    public void autoscaleBounds() {
        double[] bounds = calcAutoScale();
        if (bounds != null) {
            xLowerBoundTextField.setText(Double.toString(bounds[0]));
            xUpperBoundTextField.setText(Double.toString(bounds[1]));
            yLowerBoundTextField.setText(Double.toString(bounds[2]));
            yUpperBoundTextField.setText(Double.toString(bounds[3]));
            if (simControls instanceof CESTControls) {
                ((CESTControls) simControls).updateDeltaLimits(bounds[0], bounds[1]);
            }
            setBounds();
        }

    }

    public void autoscaleX(boolean autoX) {
        xychart.xAxis.setAutoRanging(autoX);
    }

    public void autoscaleY(boolean autoY) {
        xychart.yAxis.setAutoRanging(autoY);
    }

    public void includeXZero(boolean xZero) {
        xychart.xAxis.setZeroIncluded(xZero);
    }

    public void includeYZero(boolean yZero) {
        xychart.yAxis.setZeroIncluded(yZero);
    }

    public void updateChartEquations(String equationName, double[] pars, double[] errs, double[] fields) {
        List<GUIPlotEquation> equations = new ArrayList<>();
        String expType = getFittingMode();
        for (double field : fields) {
            double[] extras = {field / fields[0]};
            GUIPlotEquation plotEquation = new GUIPlotEquation(expType, equationName, pars, errs, extras);
            equations.add(plotEquation);
        }
        showEquations(equations);
    }

    public void showEquations(List<GUIPlotEquation> equations) {
        xychart.setEquations(equations);
//        Optional<Double> rms = rms();

    }

    public void fitR1R2NOEModel() {
        modelFitter = new FitR1R2NOEModel();
        fitIsotropicModel(modelFitter, "");
    }

    public void fitDeuteriumModel() {
        modelFitter = new FitDeuteriumModel();
        fitJCheckBox.setSelected(true);
        if (modelCheckBoxes.stream().filter(CheckBox::isSelected).findAny().isEmpty()) {
            modelCheckBoxes.stream().filter(m -> m.getText().equals("1f")).findAny().ifPresent(m -> m.setSelected(true));
        }
        fitIsotropicModel(modelFitter, "D");
    }

    public void fitIsotropicModel(FitModel fitModel, String prefix) {
        double lambdaS = sLambdaSlider.getValue();
        double lambdaTau = tauLambdaSlider.getValue();

        String tauText = tauCalcField.getText();
        Double tau = null;
        if (!tauText.isBlank()) {
            try {
                tau = Double.parseDouble(tauText);
            } catch (NumberFormatException nfE) {
            }
        } else if (prefix.equals("D")) {
            tau = 10.0;
        }
        DataIO.clearOrderPars();
        DataIO.clearSpectralDensities();
        boolean fitJ = fitJCheckBox.isSelected();
        fitModel.setLambdaS(lambdaS);
        fitModel.setLambdaTau(lambdaTau);
        fitModel.setUseLambda(lambdaCheckBox.isSelected());
        fitModel.setTau(tau);
        double tauFraction = tauFractionSlider.getValue();
        double t2Limit = t2LimitSlider.getValue();
        boolean fitTau = tauFraction > 0.001;
        fitModel.setFitTau(fitTau);
        fitModel.setT2Limit(t2Limit);
        fitModel.setNReplicates((int) nReplicatesSlider.getValue());
        fitModel.setFitJ(fitJ);
        fitModel.setBootstrapMode(bootStrapChoice.getValue());
        fitModel.setTauFraction(tauFraction);
        var modelNames = new ArrayList<String>();
        for (var modelCheckBox : modelCheckBoxes) {
            if (modelCheckBox.isSelected()) {
                String modelName = prefix + modelCheckBox.getText();
                modelNames.add(modelName);
            }
        }
        try {
            fitModel.setup(null, modelNames);
            fitModel.updaters(this::updateFitProgress, this::updateStatus);
            // fitModel.testIsoModel();
            fitModel.fitResidues();
        } catch (IllegalStateException iaE) {
            GUIUtils.warn("Model Fit Error", iaE.getMessage());
            return;
        }
    }

    public void finishModelFreeFit() {
        Double tauFit = modelFitter.getTau();
        if (tauFit != null) {
            tauCalcField.setText(String.format("%.2f", tauFit));
        }
        addMoleculeDataToAxisMenu();
        showModelFreeData();
        nReplicatesSlider.setValue(modelFitter.getNReplicates());
    }

    public void estimateCorrelationTime() {
        String r1SetName = t1Choice.getValue();
        String r2SetName = t2Choice.getValue();
        ValueSet valueSet1 = ChartUtil.getResidueProperty(r1SetName);
        ValueSet valueSet2 = ChartUtil.getResidueProperty(r2SetName);
        Map<String, Double> result = Collections.EMPTY_MAP;
        if (valueSet1 instanceof ExperimentSet) {
            if (valueSet2 instanceof ExperimentSet) {
                ExperimentSet r1Set = (ExperimentSet) valueSet1;
                ExperimentSet r2Set = (ExperimentSet) valueSet2;
                result = CorrelationTime.estimateTau(r1Set, r2Set);
            }
        } else {
            FitR1R2NOEModel fitR1R2NOEModel = new FitR1R2NOEModel();
            result = fitR1R2NOEModel.estimateTau();
        }

        if (!result.isEmpty()) {
            r1MedianField.setText(String.format("%.3f", result.get("R1")));
            r2MedianField.setText(String.format("%.3f", result.get("R2")));
            tauCalcField.setText(String.format("%.2f", result.get("tau")));
        }
    }

    public ResidueChart getActiveChart() {
        return activeChart;
    }

    public ResidueChart addChart() {
        ResidueChart newChart = ResidueChart.buildChart(barPlotCanvas);
        barCharts.add(newChart);
        activeChart = newChart;
        resizeBarPlotCanvas();
        return activeChart;
    }

    public void removeChart() {
        if ((activeChart != null) && (barCharts.size() > 1)) {
            barCharts.remove(activeChart);
            activeChart = barCharts.get(0);
            resizeBarPlotCanvas();
        }
    }

    public void updateTable(List<ExperimentData> experimentalDataSets, String fitMode) {
        ObservableList<ExperimentData.DataValue> data = FXCollections.observableArrayList();
        for (ExperimentData experimentalData : experimentalDataSets) {
            data.addAll(experimentalData.getDataValues());
        }
        resInfoTable.itemsProperty().setValue(data);

        TableColumn<ExperimentData.DataValue, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ExperimentData.DataValue, String> resColumn = new TableColumn<>("Residue");
        TableColumn<ExperimentData.DataValue, String> resNameColumn = new TableColumn<>("ResName");
        TableColumn<ExperimentData.DataValue, String> atomNameColumn = new TableColumn<>("AtomName");
        TableColumn<ExperimentData.DataValue, String> errColumn = new TableColumn<>("Error");
        TableColumn<ExperimentData.DataValue, String> peakColumn = new TableColumn<>("Peak");

        nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
        resColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
        resNameColumn.setCellValueFactory(new PropertyValueFactory<>("ResName"));
        atomNameColumn.setCellValueFactory(new PropertyValueFactory<>("AtomName"));
        errColumn.setCellValueFactory(new PropertyValueFactory<>("Error"));
        peakColumn.setCellValueFactory(new PropertyValueFactory<>("Peak"));

        resInfoTable.getColumns().clear();
        resInfoTable.getColumns().addAll(nameColumn, resNameColumn, resColumn, atomNameColumn, errColumn, peakColumn);

        if (fitMode.equals("cpmg")) {
            TableColumn<ExperimentData.DataValue, Double> xColumn = new TableColumn<>("Vcpmg");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("Reff");

            xColumn.setCellValueFactory(new PropertyValueFactory<>("X0"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resNameColumn, resColumn, atomNameColumn, xColumn, yColumn, errColumn, peakColumn);
        } else if (fitMode.equals("exp")) {
            TableColumn<ExperimentData.DataValue, Double> xColumn = new TableColumn<>("Delay");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("Intensity");

            xColumn.setCellValueFactory(new PropertyValueFactory<>("X0"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resNameColumn, resColumn, atomNameColumn,
                    //                    t1Column, t2Column, t1RhoColumn,
                    xColumn, yColumn, errColumn, peakColumn);
        } else if (fitMode.equals("cest")) {
            TableColumn<ExperimentData.DataValue, Double> x0Column = new TableColumn<>("Offset");
            TableColumn<ExperimentData.DataValue, Double> x1Column = new TableColumn<>("B1 Field");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("Intensity");

            x0Column.setCellValueFactory(new PropertyValueFactory<>("X0"));
            x1Column.setCellValueFactory(new PropertyValueFactory<>("X1"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resNameColumn, resColumn, atomNameColumn, x0Column, x1Column, yColumn, errColumn, peakColumn);
        } else if (fitMode.equals("r1rho")) {
            TableColumn<ExperimentData.DataValue, Double> x0Column = new TableColumn<>("Offset");
            TableColumn<ExperimentData.DataValue, Double> x1Column = new TableColumn<>("B1 Field");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("Intensity");

            x0Column.setCellValueFactory(new PropertyValueFactory<>("X0"));
            x1Column.setCellValueFactory(new PropertyValueFactory<>("X1"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resNameColumn, resColumn, atomNameColumn, x0Column, x1Column, yColumn, errColumn, peakColumn);
        } else if (fitMode.equals("modelfree")) {
            TableColumn<ExperimentData.DataValue, Double> xColumn = new TableColumn<>("\u03C9 (1/ns)");
            TableColumn<ExperimentData.DataValue, Double> yColumn = new TableColumn<>("log10[J(\u03C9)/1ns]");

            xColumn.setCellValueFactory(new PropertyValueFactory<>("X0"));
            yColumn.setCellValueFactory(new PropertyValueFactory<>("Y"));

            resInfoTable.getColumns().clear();
            resInfoTable.getColumns().addAll(nameColumn, resNameColumn, resColumn, atomNameColumn, xColumn, yColumn, errColumn, peakColumn);

        }
    }

    //    public void updateTableWithPars(String mapName, String[] residues, String equationName, String state, List<int[]> allStates) {
//        updateTableWithPars(mapName, residues, equationName, state, allStates, true);
//    }
//                updateTableWithPars(currentMapName, currentResidues, equationName, currentState, useStates, false);
    public void updateTableWithPars(ChartInfo chartInfo) {
        List<ParValueInterface> allParValues = new ArrayList<>();
        if (chartInfo.hasResidues() && chartInfo.hasExperiments()) {
            for (ResonanceSource resSource : chartInfo.getResidues()) {
                ExperimentResult resInfo = ChartUtil.getResInfo(chartInfo.mapName.get(0), resSource);
                if (resInfo != null) {
                    chartInfo.experimentalResult = resInfo;
                    final String useEquationName;
                    if (chartInfo.equationName.equals("best") || chartInfo.equationName.equals("+")) {
                        useEquationName = resInfo.getBestEquationName();
                    } else {
                        useEquationName = chartInfo.equationName;
                    }
                    List<ParValueInterface> parValues = resInfo.getParValues(useEquationName, chartInfo.state);
                    if (resSource.equals(chartInfo.currentResidues[0])) {
                        simControls.updateStates(chartInfo.currentStates);
                        simControls.updateSliders(parValues, useEquationName);
                        getCurrentExperimentSet().getExperimentData().stream().findFirst().ifPresent(e -> simControls.setNucleus(e.getNucleusName()));
                    }

                    allParValues.addAll(parValues);
                    CurveFit curveSet = chartInfo.experimentalResult.getCurveFit(useEquationName, chartInfo.state.replace("*", "0"));
                    if (curveSet != null) {
                        Double aic = curveSet.getParMap().get("AIC");
                        Double aicc = curveSet.getParMap().get("AICc");
                        Double rms = curveSet.getParMap().get("RMS");
                        Double rChiSq = curveSet.getParMap().get("rChiSq");
                        updateFitQuality(aic, aicc, rms, rChiSq, null);
                    }
                }
            }
            updateTableWithPars(allParValues);
        }
    }

    void updateFitQuality(Double aicValue, Double aiccValue, Double rmsValue, Double rChiSqValue, Integer nValue) {
        try {
            String aic = aicValue != null && Double.isFinite(aicValue) ? String.format("%.2f", aicValue) : "";
            String aicc = aiccValue != null && Double.isFinite(aiccValue) ? String.format("%.2f", aiccValue) : "";
            String rms = rmsValue != null && Double.isFinite(rmsValue) ? String.format("%.2f", rmsValue) : "";
            String rChiSq = rChiSqValue != null && Double.isFinite(rChiSqValue) ? String.format("%.2f", rChiSqValue) : "";
            String n = nValue != null ? String.valueOf(nValue) : "";

            aicLabel.setText(aic);
            aiccLabel.setText(aicc);
            rmsLabel.setText(rms);
            rChiSqLabel.setText(rChiSq);
            nLabel.setText(n);
        } catch (NullPointerException npEaic) {

        }

    }

    public void updateEquation(String mapName, ResonanceSource[] residues, String equationName) {
        String setEquation = "";
        for (ResonanceSource resNum : residues) {
            final String useEquationName;
            ExperimentResult resInfo = ChartUtil.getResInfo(mapName, resNum);
            if (resInfo != null) {
                if (equationName.equals("best")) {
                    useEquationName = resInfo.getBestEquationName();
                } else {
                    useEquationName = equationName;
                }

                if (setEquation.equals("")) {
                    setEquation = useEquationName;
                } else if (!setEquation.equals(useEquationName)) {
                    setEquation = "+";
                }
            }
        }
        equationChoice.setUserData("disabled");
        equationChoice.setValue(setEquation);
        equationChoice.setUserData(null);
    }

    public void updateTableWithPars(List<ParValueInterface> parValues) {
        DecimalFormat df = new DecimalFormat();
        TableColumn<ParValueInterface, String> residueColumn = new TableColumn<>("Residue");
        TableColumn<ParValueInterface, String> resNameColumn = new TableColumn<>("ResName");
        TableColumn<ParValueInterface, String> atomNameColumn = new TableColumn<>("AtomName");
        TableColumn<ParValueInterface, String> stateColumn = new TableColumn<>("State");
        TableColumn<ParValueInterface, String> nameColumn = new TableColumn<>("Name");
        TableColumn<ParValueInterface, String> valueColumn = new TableColumn<>("Value");
        TableColumn<ParValueInterface, String> errorColumn = new TableColumn<>("Error");

        residueColumn.setCellValueFactory(new PropertyValueFactory<>("Residue"));
        resNameColumn.setCellValueFactory(new PropertyValueFactory<>("ResName"));
        atomNameColumn.setCellValueFactory(new PropertyValueFactory<>("AtomName"));
        stateColumn.setCellValueFactory(new PropertyValueFactory<>("State"));
        nameColumn.setCellValueFactory(new PropertyValueFactory<>("Name"));
        valueColumn.setCellValueFactory(c -> {
            SimpleStringProperty property = new SimpleStringProperty();
            double errValue = c.getValue().getError();
            int nSig = (int) Math.floor(Math.log10(errValue)) - 1;
            nSig = -nSig;
            if (nSig < 2) {
                nSig = 2;
            }
            df.setMaximumFractionDigits(nSig);
            property.setValue(df.format(c.getValue().getValue()));
            return property;
        });
        errorColumn.setCellValueFactory(c -> {
            SimpleStringProperty property = new SimpleStringProperty();
            double errValue = c.getValue().getError();
            int nSig = (int) Math.floor(Math.log10(errValue));
            nSig = -nSig;
            if (nSig < 2) {
                nSig = 2;
            }
            df.setMaximumFractionDigits(nSig);
            property.setValue(df.format(errValue));
            return property;
        });

        parameterTable.getColumns().clear();
        parameterTable.getColumns().addAll(resNameColumn, residueColumn, atomNameColumn, stateColumn, nameColumn, valueColumn, errorColumn);
        ObservableList<ParValueInterface> data = FXCollections.observableArrayList();

        for (ParValueInterface parValue : parValues) {
            String state = parValue.getState();
            String parName = parValue.getName();
            boolean ok = true;
            // remove redundant parameters
            if ((state != null) && (parName != null)) {
                if (!state.equals("0:0:0") && (parName.equals("Kex") || parName.equals("dPPM")
                        || parName.equals("pA") || parName.equals("dPPMmin"))) {
                    ok = false;
                }
            }
            if (ok) {
                data.add(parValue);
            }
        }
        parameterTable.itemsProperty().setValue(data);
    }

    public void selectTableRow(String seriesName, int index) {
        parTabPane.getSelectionModel().select(1);
        List<ExperimentData.DataValue> data = resInfoTable.getItems();
        int iRow = 0;
        for (ExperimentData.DataValue dValue : data) {
            if ((dValue.getIndex() == index) && ((dValue.getName() + ":" + dValue.getResidue()).equals(seriesName))) {
                resInfoTable.getSelectionModel().clearAndSelect(iRow);
                resInfoTable.scrollTo(iRow);
                return;
            }
            iRow++;

        }

    }

    public String getPeakNumFromTable() { //getPeakNumFromTable(String seriesName, int index)
        parTabPane.getSelectionModel().select(1);
        List<ExperimentData.DataValue> data = resInfoTable.getItems();
        int iRow = 0;
        int peakNum;
        String name;
        name = data.get(iRow).getName();
        peakNum = data.get(iRow).getPeak();
        return name + "." + peakNum;
    }

    public void clearChart() {
        activeChart.getData().clear();
    }

    public void showRelaxationValues(String setName, String valueName, String parName) {
        ValueSet valueSet = ChartUtil.getResidueProperty(setName);
        if (valueSet != null) {
            List<RelaxationValues> values = (List<RelaxationValues>) valueSet.rValues();
            ObservableList<DataSeries> data = ChartUtil.getRelaxationDataSeries(values, valueName, setName, parName);
            String yLabel = valueName.equalsIgnoreCase(parName) ? parName
                    : valueName.toUpperCase() + ": " + parName;

            addSeries(data, setName, yLabel, true);
        }

    }

    public void setYAxisType(String expMode, String setName, String eqnName, String state, String typeName, boolean updateProps) {
        ObservableList<DataSeries> data = ChartUtil.getParMapData(setName, eqnName, state, typeName);
        String yLabel = expMode.equalsIgnoreCase(typeName) ? typeName
                : expMode.toUpperCase() + ": " + typeName;
        addSeries(data, setName, yLabel, updateProps);
    }

    public void addSeries(ObservableList<DataSeries> data, String setName, String yLabel, boolean updateProps) {
        ResidueChart chart = activeChart;
        if (updateProps) {
            chart.setResProps(ChartUtil.getResidueProperty(setName));
        }

        chart.getData().addAll(data);
        int i = 0;
        for (DataSeries series : chart.getData()) {
            series.setFill(PlotData.colors[i++ % PlotData.colors.length]);
        }
        chart.currentSeriesName = data.get(0).getName();
        Axis xAxis = chart.xAxis;
        Axis yAxis = chart.yAxis;
        xAxis.setAutoRanging(false);
        chart.autoScale(false);
        double xMin = Math.floor((ChartUtil.minRes - 2) / 5.0) * 5.0;
        double xMax = Math.ceil((ChartUtil.maxRes + 2) / 5.0) * 5.0;
        barChartXMin = Optional.of(xMin);
        barChartXMax = Optional.of(xMax);
        xAxis.setLowerBound(xMin);
        xAxis.setUpperBound(xMax);
        if (yLabel.startsWith("S") && yLabel.endsWith("2")) {
            yAxis.setAutoRanging(false);
            yAxis.setLowerBound(0.0);
            yAxis.setUpperBound(1.0);
        } else {
            yAxis.setAutoRanging(true);
        }
        chart.yAxis.setLabel(yLabel);

        chart.setTitle(setName);
        setCurrentExperimentSet(ChartUtil.getResidueProperty(setName));
        refreshResidueCharts();
    }

    void makeT1T2Menu() {
        t1Choice.getItems().clear();
        t2Choice.getItems().clear();
        // noeChoice.getItems().clear();
        t1Choice.getItems().add("");
        t2Choice.getItems().add("");
        //noeChoice.getItems().add("");
        Collection<String> setNames = ChartUtil.getResiduePropertyNames();
        for (var setName : setNames) {
            ValueSet valueSet = ChartUtil.getResidueProperty(setName);
            if (valueSet instanceof ExperimentSet) {
                t1Choice.getItems().add(setName);
                t2Choice.getItems().add(setName);
            }
        }
    }

    void makeRelaxSetMenu() {
        relaxSetMenu.getItems().clear();
        MoleculeBase moleculeBase = MoleculeFactory.getActive();
        Set<Integer> fieldSets = new HashSet<>();
        if (moleculeBase != null) {
            var relaxSets = moleculeBase.relaxationSetMap();
            for (var entry : relaxSets.entrySet()) {
                int field = (int) entry.getValue().field();
                if (!fieldSets.contains(field)) {
                    CheckMenuItem checkMenuItem = new CheckMenuItem(String.valueOf(field));
                    checkMenuItem.setSelected(entry.getValue().active());
                    checkMenuItem.setOnAction(e -> activateRelaxField(field, checkMenuItem.isSelected()));
                    relaxSetMenu.getItems().add(checkMenuItem);
                    fieldSets.add(field);
                }
            }
        }
    }

    void activateRelaxField(int field, boolean state) {
        MoleculeBase moleculeBase = MoleculeFactory.getActive();
        if (moleculeBase != null) {
            var relaxSets = moleculeBase.relaxationSetMap();
            for (var entry : relaxSets.entrySet()) {
                int entryField = (int) entry.getValue().field();
                if (entryField == field) {
                    entry.getValue().active(state);
                }
            }
        }
    }

    void makeAxisMenu() {
        makeT1T2Menu();
        experimentalDataAxisMenu.getItems().clear();
        moleculeDataAxisMenu.getItems().clear();
        addMoleculeDataToAxisMenu();
        addResiduePropertiesToAxisMenu();
    }

    void addMoleculeDataToAxisMenu() {
        addRelaxationDataToAxisMenu(DataIO.getRelaxationDataFromMolecule());
        addOrderParDataToAxisMenu(DataIO.getOrderParSetFromMolecule());
    }

    void addRelaxationDataToAxisMenu(Map<String, RelaxationSet> molResProps) {
        for (var entry : molResProps.entrySet()) {
            String setName = entry.getKey();
            ChartUtil.addResidueProperty(setName, molResProps.get(setName));
            RelaxationSet relaxSet = entry.getValue();
            Menu cascade = new Menu(setName);
            var values = relaxSet.values();
            if (!values.isEmpty()) {
                RelaxationValues value = values.get(0);
                moleculeDataAxisMenu.getItems().add(cascade);
                String[] parNames = value.getParNames();
                for (var parName : parNames) {
                    MenuItem cmItem1 = new MenuItem(parName);
                    cmItem1.setOnAction(e -> showRelaxationValues(setName, setName, parName));
                    cascade.getItems().add(cmItem1);
                }

            }
        }
    }

    void addOrderParDataToAxisMenu(Map<String, OrderParSet> molResProps) {
        for (var entry : molResProps.entrySet()) {
            String setName = entry.getKey();
            ChartUtil.addResidueProperty(setName, molResProps.get(setName));
            OrderParSet relaxSet = entry.getValue();
            Menu cascade = new Menu(setName);
            var values = relaxSet.values();
            if (!values.isEmpty()) {
                RelaxationValues value = values.get(0);
                moleculeDataAxisMenu.getItems().add(cascade);
                String[] parNames = value.getParNames();
                for (var parName : parNames) {
                    MenuItem cmItem1 = new MenuItem(parName);
                    cmItem1.setOnAction(e -> showRelaxationValues(setName, setName, parName));
                    cascade.getItems().add(cmItem1);
                }
            }
        }
    }

    void addOrderParSetsToAxisMenu() {
        System.out.println("add order");
        Map<String, OrderParSet> molResProps = DataIO.getOrderParSetFromMolecule();
        orderParSetAxisMenu.getItems().clear();
        for (var entry : molResProps.entrySet()) {
            String setName = entry.getKey();
            if (!entry.getValue().values().isEmpty()) {
                CheckMenuItem checkMenuItem = new CheckMenuItem(setName);
                checkMenuItem.setSelected(entry.getValue().active());
                checkMenuItem.setOnAction(e -> activateOrderParSet(entry.getValue(), checkMenuItem.isSelected()));
                orderParSetAxisMenu.getItems().add(checkMenuItem);
            }
        }
    }

    private void activateOrderParSet(OrderParSet orderParSet, boolean state) {
        orderParSet.active(state);
    }

    public Map<String, ResidueChart> setupCharts(List<String> chartNames) {
        int nNewCharts = chartNames.size();
        int nCharts = barCharts.size();
        if (nCharts > nNewCharts) {
            barCharts.subList(nNewCharts, nCharts).clear();
        }
        nCharts = barCharts.size();
        for (int iChart = nCharts; iChart < nNewCharts; iChart++) {
            addChart();
        }
        int iChart = 0;
        var chartMap = new HashMap<String, ResidueChart>();
        for (ResidueChart residueChart : barCharts) {
            activeChart = residueChart;
            clearChart();
            chartMap.put(chartNames.get(iChart), residueChart);
            iChart++;
            boolean xVisible = iChart == barCharts.size();
            residueChart.getXAxis().setLabelVisible(xVisible);
            residueChart.getXAxis().setTickLabelsVisible(xVisible);
        }
        return chartMap;

    }

    public void showT1T2NoeData() {
        var molResProps = DataIO.getRelaxationDataFromMolecule();

        List<String> chartNames = molResProps.values().stream().
                filter(v -> !v.values().isEmpty()).
                filter(v -> v.values().get(0) != null).
                filter(v -> v.active()).
                map(v -> v.relaxType().getName()).
                collect(Collectors.toSet()).stream().sorted().collect(Collectors.toList());

        var chartMap = setupCharts(chartNames);
        molResProps.values().stream().
                filter(v -> !v.values().isEmpty()).
                filter(v -> v.values().get(0) != null).
                filter(v -> v.active()).
                sorted((a, b) -> {
                    double fieldA = a.field();
                    double fieldB = b.field();
                    return Double.compare(fieldA, fieldB);
                }).
                forEach(v -> {
                    var setName = v.name();
                    var rData = (RelaxationData) v.values().get(0);
                    String[] parNames = rData.getParNames();
                    activeChart = chartMap.get(v.relaxType().getName());
                    showRelaxationValues(setName, parNames[0], parNames[0]);
                });
        sortChartSeries();
        resizeBarPlotCanvas();
    }

    public void showModelFreeData() {
        var chartNames = List.of("S2", "Sf2", "Ss2", "Tau_e",
                "Tau_f", "Tau_s", "Rex", "model", "rchisq");
        showModelFreeData(chartNames);
    }

    public void showModelFreeDataS2() {
        var chartNames = List.of("S2", "Sf2", "Ss2");
        showModelFreeData(chartNames);
    }

    public void showModelFreeDataTau() {
        var chartNames = List.of("Tau_e", "Tau_f", "Tau_s");
        showModelFreeData(chartNames);
    }

    public void showModelFreeData(List<String> chartNames) {
        var chartMap = setupCharts(chartNames);
        var usedSet = new TreeSet<String>();
        var molResProps = DataIO.getOrderParSetFromMolecule();
        molResProps.entrySet().stream().filter(entry -> entry.getValue().active()).
                forEach(entry -> {
                    OrderParSet v = entry.getValue();
                    String key = entry.getKey();
                    var values = v.values();
                    boolean hasNull = values.stream().anyMatch(value -> value.getValue() == null);
                    if (!hasNull) {
                        var setName = v.name();
                        for (var parName : chartNames) {
                            boolean hasValue = values.stream().anyMatch(value
                                    -> (value.getValue(parName) != null) && (value.getValue(parName) > 1.0e-6));
                            if (hasValue) {
                                activeChart = chartMap.get(parName);
                                showRelaxationValues(setName, parName, parName);
                                usedSet.add(parName);
                            }
                        }
                    }
                });
        for (String chartName : chartNames) {
            if (!usedSet.contains(chartName)) {
                ResidueChart resChart = chartMap.get(chartName);
                barCharts.remove(resChart);
            }
        }
        sortChartSeries();
        resizeBarPlotCanvas();
    }

    private void sortChartSeries() {
        for (ResidueChart residueChart : barCharts) {
            Collections.sort(residueChart.getData(), seriesComparator);
        }
    }

    class SeriesComparator implements Comparator<DataSeries> {

        @Override
        public int compare(DataSeries s1, DataSeries s2) {
            String o1 = s1.getName();
            String o2 = s2.getName();
            String[] fields1 = o1.split("\\|");
            String[] fields2 = o2.split("\\|");
            String a1 = fields1[fields1.length - 1];
            String a2 = fields2[fields2.length - 1];
            int c = a1.compareTo(a2);
            if (c == 0) {
                c = o1.compareTo(o2);
            }
            return c;
        }
    }

    void addResiduePropertiesToAxisMenu() {
        Collection<String> setNames = ChartUtil.getResiduePropertyNames();
        setNames.stream().sorted().forEach(setName -> {
            var valueSet = ChartUtil.getResidueProperty(setName);
            if (valueSet instanceof ExperimentSet) {
                ExperimentSet experimentSet = (ExperimentSet) valueSet;
                Menu cascade = new Menu(setName);
                experimentalDataAxisMenu.getItems().add(cascade);
                String expMode = experimentSet.getExpMode();
                String[] parTypes = getParTypes(experimentSet.getExpMode());
                if (experimentSet.getEquationNames().size() == 0) {
                    final String parName;
                    if (expMode.equals("r1") || expMode.equals("r2") || expMode.equals("rq") || expMode.equals("rap")) {
                        parName = "R";
                    } else if (experimentSet.getExpMode().equals("noe")) {
                        parName = "NOE";
                    } else {
                        parName = "Kex";
                    }
                    MenuItem cmItem1 = new MenuItem("R");
                    cmItem1.setOnAction(e -> setYAxisType(expMode, experimentSet.name(), "best", "0:0:0", parName, true));
                    cascade.getItems().add(cmItem1);
                } else {
                    for (String parType : parTypes) {
                        if (experimentSet.getEquationNames().size() == 1) {
                            String equationName = experimentSet.getEquationNames().get(0);
                            MenuItem cmItem1 = new MenuItem(parType);
                            cmItem1.setOnAction(e -> setYAxisType(expMode, setName, equationName, "0:0:0", parType, true));
                            cascade.getItems().add(cmItem1);
                        } else {

                            Menu cascade2 = new Menu(parType);
                            cascade.getItems().add(cascade2);
                            ArrayList<String> equationNames = new ArrayList<>();
                            if (experimentSet.getEquationNames().size() > 1) {
                                equationNames.add("best");
                            }
                            equationNames.addAll(experimentSet.getEquationNames());
                            List<String> stateStrings = experimentSet.getStateStrings();
                            if (equationNames.size() == 0) {
                                MenuItem cmItem1 = new MenuItem(experimentSet.name());
                                cmItem1.setOnAction(e -> setYAxisType(expMode, experimentSet.name(), "best", "0:0:0", "R", true));
                                cascade2.getItems().add(cmItem1);

                            } else {

                                for (String equationName : equationNames) {
                                    if ((stateStrings.size() < 2) || parType.equals("RMS") || parType.equals("AIC") || parType.equals("Equation")) {
                                        MenuItem cmItem1 = new MenuItem(equationName);
                                        cmItem1.setOnAction(e -> setYAxisType(expMode, setName, equationName, "0:0:0", parType, true));
                                        cascade2.getItems().add(cmItem1);
                                    } else {
                                        boolean isValidPar = false;
                                        if (equationName.equals("best")) {
                                            isValidPar = true;
                                        } else {
                                            String[] validPars = getParNames(equationName);
                                            for (String validPar : validPars) {
                                                if (validPar.equals(parType)) {
                                                    isValidPar = true;
                                                    break;
                                                }
                                            }
                                        }
                                        if (isValidPar) {
                                            Menu cascade3 = new Menu(equationName);
                                            cascade2.getItems().add(cascade3);
                                            experimentSet.getStateStrings().forEach(state -> {
                                                MenuItem cmItem1 = new MenuItem(state);
                                                cmItem1.setOnAction(e -> setYAxisType(expMode, setName, equationName, state, parType, true));
                                                cascade3.getItems().add(cmItem1);

                                            });
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        });

    }

    @FXML
    public void guesses() {
        if (!hasExperimentSet()) {
            guessSimData();
        } else {
            try {
                EquationFitter equationFitter = getFitter();
                ResonanceSource[] resNums = {chartInfo.getSource()};
                equationFitter.setData(getCurrentExperimentSet(), resNums);
                String equationName = simControls.getEquation();
                List<ParValueInterface> guesses = equationFitter.guessPars(equationName);
                if (guesses != null) {
                    simControls.updateSliders(guesses, equationName);
                    simControls.simSliderAction("");
                }
            } catch (NullPointerException npE1) {
                Alert alert = new Alert(Alert.AlertType.ERROR);
                alert.setContentText("Error: Residue must be selected");
                alert.showAndWait();
            }
        }
    }

    public Optional<Double> rms() {
        Optional<Double> rms = Optional.empty();
        if (hasExperimentSet()) {
            try {
                EquationFitter equationFitter = getFitter();
                if (chartInfo.hasResult()) {
                    ResonanceSource[] resNums = {chartInfo.getSource()};
                    equationFitter.setData(getCurrentExperimentSet(), resNums);
                    String equationName = simControls.getEquation();
                    equationFitter.setupFit(equationName);
                    int[][] map = equationFitter.getFitModel().getMap();
                    double[] sliderGuesses = simControls.sliderGuess(equationName, map);
                    rms = Optional.of(equationFitter.rms(sliderGuesses));
                } else {
                    rms = Optional.empty();
                }
            } catch (NullPointerException npE2) {
                rms = Optional.empty();
            }
        }
        return rms;
    }

    @FXML
    public void fitEquation() {
        fitResult = null;
        try {
            EquationFitter equationFitter = getFitter();
            if (!hasExperimentSet()) {
                fitSimData();
            } else {
                ResonanceSource[] resNums = {chartInfo.getSource()};
                equationFitter.setData(getCurrentExperimentSet(), resNums);
                String equationName = simControls.getEquation();
                equationFitter.setupFit(equationName);
                int[][] map = equationFitter.getFitModel().getMap();
                double[] sliderGuesses = null;
                if (sliderGuessCheckBox.isSelected()) {
                    sliderGuesses = simControls.sliderGuess(equationName, map);
                }
                CoMDOptions options = new CoMDOptions(true);
                fitResult = equationFitter.doFit(equationName, sliderGuesses, options);
                updateAfterFit(fitResult);
            }
        } catch (NullPointerException npE2) {
            npE2.printStackTrace();
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setContentText("Error: Residue must be selected when fitting");
            alert.showAndWait();
        }
    }

    public void updateAfterFit(FitResult fitResult) {
        List<GUIPlotEquation> equations = new ArrayList<>();
        int nCurves = fitResult.getNCurves();
        for (int iCurve = 0; iCurve < nCurves; iCurve++) {
            CurveFit curveFit = fitResult.getCurveFit(iCurve);
            List<ParValueInterface> parValues = curveFit.getParValues();
            equationChoice.getSelectionModel().select(fitResult.getEquationName());
            String aic = String.format("%.2f", fitResult.getAicc());
            String rms = String.format("%.3f", fitResult.getRms());
            String rChiSq = String.format("%.2f", fitResult.getRChiSq());
            updateFitQuality(null, fitResult.getAicc(), fitResult.getRms(), fitResult.getRChiSq(), null);
            aiccLabel.setText(aic);
            rmsLabel.setText(rms);
            rChiSqLabel.setText(rChiSq);
            updateTableWithPars(parValues);
            simControls.updateSliders(parValues, fitResult.getEquationName());
            String equationName = fitResult.getEquationName(); //equationSelector.getValue();
            String expType = getFittingMode();
            Optional<Experiment> optionalData = Optional.empty();
            if (hasExperimentSet()) {
                optionalData = getCurrentExperimentSet().getExperimentData().stream().findFirst();
            }

            if (optionalData.isPresent() && optionalData.get().getExtras().size() > 0) {
                for (Experiment expData : getCurrentExperimentSet().getExperimentData()) {
                    double[] pars = curveFit.getEquation().getPars(); //pars = getPars(equationName);
                    double[] errs = curveFit.getEquation().getErrs(); //double[] errs = new double[pars.length];
                    double[] extras = new double[3];
                    extras[0] = expData.getNucleusField();
                    extras[1] = expData.getExtras().get(0);
                    extras[2] = expData.getExtras().get(1);
                    GUIPlotEquation plotEquation = new GUIPlotEquation(expType, equationName, pars, errs, extras);
                    equations.add(plotEquation);
                }
            } else {
                GUIPlotEquation plotEquation = new GUIPlotEquation(expType, curveFit.getEquation());
                equations.add(plotEquation);

            }
        }
        showEquations(equations);
    }

    @FXML
    public void fitResidues() {
        fitResult = null;
        if (hasExperimentSet()) {
            residueFitter.fitResidues(getCurrentExperimentSet());
        }
    }

    public void fitResiduesNow() {
        fitResult = null;
        if (hasExperimentSet()) {
            residueFitter.fitResiduesNow(getCurrentExperimentSet());
        }
    }

    @FXML
    public void fitGroupResidues() {
        if (hasExperimentSet()) {
            fitResult = null;
            List<List<ResonanceSource>> allResidues = new ArrayList<>();
            fittingResidues.clear();
            fittingResidues.addAll(ResidueChart.selectedResidues);
            List<ResonanceSource> groupResidues = new ArrayList<>(ResidueChart.selectedResidues);
            if (!groupResidues.isEmpty()) {
                allResidues.add(groupResidues);
                residueFitter.fitResidues(getCurrentExperimentSet(), allResidues);
            }
        }
//        }
    }

    public void refreshFit() {
        ResidueChart chart = getActiveChart();
        ResidueChart.selectedResidues.clear();
        ResidueChart.selectedResidues.addAll(fittingResidues);
        refreshResidueCharts();
        chart.showInfo();
        makeAxisMenu();
    }

    @FXML
    public void haltFit() {
        if (modelFitter != null) {
            modelFitter.haltFit();
        } else {
            residueFitter.haltFit();
            makeAxisMenu();
        }
    }

    @FXML
    public void saveParameters() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Save Parameter File");
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            DataIO.saveResultsFile(file.getAbsolutePath(), getCurrentExperimentSet(), false);
        }
    }

    @FXML
    public void saveR1R2NOE() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Save R1/R2/... File");
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            try {
                RelaxationData.writeToFile(file);
            } catch (IOException e) {
                ExceptionDialog exceptionDialog = new ExceptionDialog(e);
                exceptionDialog.showAndWait();
            }
        }
    }

    @FXML
    public void saveOrderParameters() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Save Order Parameters File");
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            try {
                OrderPar.writeToFile(file, ",");
            } catch (IOException e) {
                ExceptionDialog exceptionDialog = new ExceptionDialog(e);
                exceptionDialog.showAndWait();
            }
        }
    }

    public void saveSpectralDensities() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Save Spectral Densities File");
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            try {
                SpectralDensity.writeToFile(file, ",");
            } catch (IOException e) {
                ExceptionDialog exceptionDialog = new ExceptionDialog(e);
                exceptionDialog.showAndWait();
            }
        }
    }

    public void saveParametersSTAR() throws IOException, InvalidMoleculeException, ParseException, InvalidPeakException {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Save STAR File");
        File file = fileChooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            DataIO.writeSTAR3File(file.getAbsolutePath());
            System.out.println("wrote " + file.getAbsolutePath());
        }
    }

    public void addRelaxResultsToMol() {
        MoleculeBase mol = MoleculeFactory.getActive();
        String alertText = "Add Relax results to map?";
        if (mol != null) {
            alertText = "Add relax results to molecule " + mol.getName() + "?";
        }
        Alert alert = new Alert(Alert.AlertType.CONFIRMATION, alertText);
        Optional<ButtonType> response = alert.showAndWait();
        if (response.isPresent() && response.get().getText().equals("OK")) {
            Collection<String> setNames = ChartUtil.getResiduePropertyNames();
            setNames.stream().sorted().forEach(setName -> {
                var valueSet = ChartUtil.getResidueProperty(setName);
                if (valueSet instanceof ExperimentSet) {
                    ExperimentSet experimentSet = (ExperimentSet) valueSet;
                    RelaxTypes expMode = RelaxTypes.valueOf(experimentSet.getExpMode().toUpperCase());
                    DataIO.addRelaxationFitResults(experimentSet, expMode);
                }
            });
        }
    }

    public Double updateFitProgress(Double f) {
        if (hasExperimentSet()) {
            ResidueChart chart = getActiveChart();
            String seriesName = chart.currentSeriesName;
            String[] sParts;
            if (seriesName.length() == 0) {
                sParts = new String[4];
                sParts[0] = getCurrentExperimentSet().name();
                sParts[1] = "best";
                sParts[2] = "0:0:0";
                sParts[3] = "RMS";
            } else {
                sParts = seriesName.split("\\|");
            }
            if (Platform.isFxApplicationThread()) {
                clearChart();
                statusBar.setProgress(f);
                setYAxisType(getCurrentExperimentSet().getExpMode(), sParts[0], sParts[1], sParts[2], sParts[3], false);
                setCurrentExperimentSet(ChartUtil.getResidueProperty(getCurrentExperimentSet().name()));

            } else {
                Platform.runLater(() -> {
                    clearChart();
                    setYAxisType(getCurrentExperimentSet().getExpMode(), sParts[0], sParts[1], sParts[2], sParts[3], false);
                    statusBar.setProgress(f);
                    setCurrentExperimentSet(ChartUtil.getResidueProperty(getCurrentExperimentSet().name()));
                });
            }
        } else if (modelFitter != null) {
            statusBar.setProgress(f);
        }
        return null;

    }

    public void processingDone() {
        makeAxisMenu();
    }

    public Double updateStatus(ProcessingStatus status) {
        if (Platform.isFxApplicationThread()) {
            updateStatusNow(status);
        } else {
            Platform.runLater(() -> updateStatusNow(status));

        }
        return null;
    }

    public void updateStatusNow(ProcessingStatus status) {
        String s = status.getStatus();
        if (s == null) {
            statusBar.setText("");
        } else {
            statusBar.setText(s);
            if (s.equals("Done")) {
                if (modelFitter != null) {
                    finishModelFreeFit();
                } else {
                    refreshFit();
                }
            }
        }
        if (status.isOk()) {
            statusCircle.setFill(Color.GREEN);
        } else {
            if (status.getThrowable() != null) {
                status.getThrowable().printStackTrace();
            }
            statusCircle.setFill(Color.RED);
        }
        statusBar.setProgress(0.0);
    }

    public void saveBarChart() throws IOException {
        FileChooser chooser = new FileChooser();
        chooser.setInitialFileName("barchart.png");
        File file = chooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            snapit(chartBox, file);
        }
    }

    public void saveXYChart() throws IOException {
        FileChooser chooser = new FileChooser();
        chooser.setInitialFileName("cpmgchart.png");
        File file = chooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            snapit(chartPane, file);
        }
    }

    public void exportExecutable(String fileSuggestion, String exportType, boolean saveBar) throws ScriptException {
        FileChooser chooser = new FileChooser();
        chooser.setInitialFileName(fileSuggestion);
        File file = chooser.showSaveDialog(MainApp.primaryStage);
        if (file != null) {
            String filePath = file.getAbsolutePath();

            Map<String, Object>[] plottedData = xychart.getPlottedData();
            Map<String, Object> graphData;
            graphData = xychart.getGraphData();
            graphData.put("file", filePath);
            graphData.put("exportType", exportType);
            ArrayList<Object> barChartData;
            if (!"grace".equals(exportType) && saveBar) {
                barChartData = new ArrayList<>(1);
            } else {
                barChartData = new ArrayList<>(0);
            }
            var interpreter = MainApp.getInterpreter();

            interpreter.exec("from exporter import Writer");
            interpreter.set("data", plottedData);
            interpreter.set("configData", graphData);
            interpreter.set("barChartData", barChartData);
            interpreter.exec("writer = Writer(configData,data,barChartData)");
            interpreter.exec("writer.writeFromExportData()");
        }
    }

    public void saveGraceFile() throws ScriptException {
        exportExecutable("ASCII.agr", "grace", false);
    }

    public void savePythonFile() throws ScriptException {
        exportExecutable("graph.py", "python", false);
    }

    public void saveRFile() throws ScriptException {
        exportExecutable("graph.r", "r", false);
    }

    public void saveBarToPythonFile() throws ScriptException {
        exportExecutable("graph.py", "python", true);
    }

    public void saveBarToRFile() throws ScriptException {
        exportExecutable("graph.r", "r", true);
    }

    public void snapit(Node node, File file) throws IOException {
        double scale = 4.0;
        final Bounds bounds = node.getLayoutBounds();
        final WritableImage image = new WritableImage(
                (int) Math.round(bounds.getWidth() * scale),
                (int) Math.round(bounds.getHeight() * scale));
        final SnapshotParameters spa = new SnapshotParameters();
        spa.setTransform(javafx.scene.transform.Transform.scale(scale, scale));
        node.snapshot(spa, image);
        ImageIO.write(SwingFXUtils.fromFXImage(image, null), "png", file);
    }

    public FitResult getFitResult() {
        return fitResult;
    }

    public ExperimentResult getResidueInfo() {
        return chartInfo.getResult();
    }

    public String getEquationName() {
        return chartInfo.getEquationName();
    }

    public String getFittingMode() {
        String fitMode = "cpmg";
        if (!hasExperimentSet()) {
            String simMode = getSimMode();
            if (simMode != null) {
                fitMode = simMode;
            } else {
                fitMode = "modelfree";
            }
        } else {
            fitMode = getCurrentExperimentSet().getExpMode().toLowerCase();
            if (fitMode.equals("r1") || fitMode.equals("r2") || fitMode.equals("rq") || fitMode.equals("rap")) {
                fitMode = "exp";
            }
        }
        return fitMode.toLowerCase();
    }

    public EquationFitter getFitter() {
        CoMDOptions options = new CoMDOptions(true);
        if (getFittingMode().equals("exp")) {
            return new ExpFitter(options);
        } else if (getFittingMode().equals("cpmg")) {
            return new CPMGFitter(options);
        } else if (getFittingMode().equals("cest")) {
            return new CESTFitter(options);
        } else if (getFittingMode().equals("r1rho")) {
            return new R1RhoFitter(options);
        }
        return null;
    }

    public String[] getParNames(String equationName) {
        String fitMode = getFittingMode();
        EquationType type = ResidueFitter.getEquationType(fitMode, equationName);
        return type.getParNames();
    }

    public String[] getParTypes(String mode) {
        String[] cpmgTypes = {"R2", "dPPMmin", "Kex", "pA", "dPPM", "RMS", "AIC", "Equation"};
        String[] expTypes = {"A", "R", "C", "RMS", "AIC", "Equation"};
        String[] sTypes = {"S2", "Rex"};
        String[] cestTypes = {"kex", "pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B", "RMS", "AIC", "Equation"};
        String[] r1rhoTypes = {"kex", "pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B", "RMS", "AIC", "Equation"};
        String[] nullTypes = {"RMS", "AIC", "Equation"};
        String[] noeTypes = {"NOE"};
        switch (mode) {
            case "exp":
                return expTypes;
            case "cpmg":
                return cpmgTypes;
            case "cest":
                return cestTypes;
            case "r1rho":
                return r1rhoTypes;
        }
        switch (mode) {
            case "r1":
            case "r2":
            case "rq":
            case "rap":
                return expTypes;
            case "s2":
                return sTypes;
            case "noe":
                return noeTypes;
        }
        return nullTypes;
    }

    @FXML
    void setBestEquation() {
        for (ResonanceSource resSource : chartInfo.currentResidues) {
            ExperimentResult resInfo = ChartUtil.getResInfo(chartInfo.mapName.get(0), resSource);
            if (resInfo != null) {
                String equationName = equationChoice.getValue();
                resInfo.setBestEquationName(equationName);
            }
        }
        updateFitProgress(1.0);

    }

    void equationAction() {
        System.out.println("eq act " + equationChoice.getUserData());
        if (equationChoice.getUserData() == null) {
            String equationName = equationChoice.getValue();
            if (!chartInfo.currentStates.isEmpty() && equationName != null) {
                // copy it so it doesn't get cleared by clear call in updateTableWithPars
                chartInfo.equationName = equationName;
                updateTableWithPars(chartInfo);
                showInfo(equationName);
            }
        }
    }

    public String getParametersEquation() {
        return equationChoice.getValue();
    }

    public String getSimMode() {
        String simSelected = simChoice.getSelectionModel().getSelectedItem();
        if (simSelected == null) {
            return "cpmg";
        }
        return simSelected.substring(9).toLowerCase();
    }

    void simAction() {
        clearProject(false);
        getSimMode();
        setSimControls();
        updateXYChartLabels();
        genDataSDevTextField.setText("");
        simControls.simSliderAction("");
    }

    void bootStrapChanged(FitModel.BootstrapMode mode) {
        if (mode != FitModel.BootstrapMode.PARAMETRIC) {
            fitJCheckBox.setSelected(true);
            int nReplicates = (int) nReplicatesSlider.getValue();
            if (nReplicates < 10) {
                nReplicates = 50;
                nReplicatesSlider.setValue(nReplicates);
            }
            if (!lambdaCheckBox.isSelected()) {
                modelCheckBoxes.forEach(checkBox -> checkBox.setSelected(true));
            }
        }
    }

    void lambdaChanged(boolean state) {
        if (state) {
            fitJCheckBox.setSelected(true);
            if (sLambdaSlider.getValue() < 1.0e-5) {
                sLambdaSlider.setValue(0.5);
                tauLambdaSlider.setValue(0.5);
            }
            modelCheckBoxes.forEach(checkBox -> checkBox.setSelected(false));
            modelCheckBoxes.get(modelCheckBoxes.size() - 1).setSelected(true);
        }
    }

    @FXML
    void clearProject() {
        if (GUIUtils.affirm("Clear Project")) {
            clearProject(true);
        }
    }

    void clearProject(boolean clearXY) {
        chartInfo.clear();
        setCurrentExperimentSet(null);
        fitResult = null;
        ChartUtil.clearResidueProperties();
        if (clearXY) {
            xychart.clear();
        } else {
            xychart.getData().clear();
        }
        clearSecondaryStructure();
        barCharts.remove(activeChart);
        addChart();
    }

    @FXML
    void clearRelaxPars() {
        if (GUIUtils.affirm("Clear Relaxation Parameters")) {
            DataIO.clearRelaxationData();
        }
    }

    @FXML
    void clearOrderPars() {
        if (GUIUtils.affirm("Clear Order Parameters")) {
            DataIO.clearOrderPars();
        }
    }

    void showInfo(String equationName) {
        if (hasExperimentSet()) {
            showInfo(chartInfo, xychart);
        }
    }

    void showInfo(ChartInfo chartInfo, PlotData plotData) {
        ArrayList<GUIPlotEquation> equations = new ArrayList<>();
        ObservableList<DataSeries> allData = FXCollections.observableArrayList();
        List<ExperimentData> experimentalDataSets = new ArrayList<>();
        List<int[]> allStates = new ArrayList<>();
        boolean calcScale = scalePlot.isSelected();
        List<ParValueInterface> parValues = null;
        String fitMode = getFittingMode();
        int iSeries = 0;
        if (chartInfo.hasExperiments() && chartInfo.hasResidues()) {
            for (Experiment expData : ((ExperimentSet) chartInfo.getExperiments()).getExperimentData()) {
                if (!ExperimentSet.matchStateString(chartInfo.state, expData.getState())) {
                    // continue;
                }
                String expName = expData.getName();
                for (ResonanceSource resNum : chartInfo.getResidues()) {
                    if (expData.getResidueData(resNum) != null) {
                        experimentalDataSets.add(expData.getResidueData(resNum));
                        DataSeries series = ChartUtil.getMapData(chartInfo.mapName.get(0), expName, resNum);
                        series.setStroke(PlotData.colors[iSeries % 8]);
                        series.setFill(PlotData.colors[iSeries % 8]);
                        allData.add(series);
                        GUIPlotEquation equation = ChartUtil.getEquation(expData,
                                chartInfo.mapName.get(0), resNum, chartInfo.equationName, expData.getState(),
                                expData.getNucleusField());
                        double maxY = 1.0;
                        if (equation != null) {
                            equation.setColor(PlotData.colors[iSeries % 8]);
                            if (calcScale) {
                                maxY = equation.calculate(equation.getMinX()) / 100.0;
                            }
                            equation.setScaleValue(maxY);
                            equations.add(equation);
                        } else if (calcScale) {
                            var maxOpt = series.getValues().stream().mapToDouble(XYValue::getYValue).max();
                            if (maxOpt.isPresent()) {
                                maxY = maxOpt.getAsDouble();
                            }
                        }
                        series.setScale(maxY);
                        iSeries++;
                    }
                }
                ExperimentSet expSet = (ExperimentSet) chartInfo.getExperiments();
                int[] states = expSet.getStateIndices(0, expData);
                allStates.add(states);
            }
        } else {
            fitMode = "modelfree";
            for (ResonanceSource resonanceSource : chartInfo.getResidues()) {
                Atom atom = resonanceSource.getAtom();
                System.out.println("ressourc " + resonanceSource + " " + atom.getFullName());
                Map<String, SpectralDensity> spectralDensityMap = atom.getSpectralDensity();
                var sdData = ChartUtil.getSpectralDensityData(spectralDensityMap);
                allData.addAll(sdData);
                List<Double> xV = new ArrayList<>();
                List<Double> yV = new ArrayList<>();
                List<Double> eV = new ArrayList<>();
                for (var v : sdData) {
                    var d = v.getData();
                    for (var dv : d) {
                        if (dv instanceof XYEValue xye) {
                            xV.add(xye.getXValue());
                            yV.add(xye.getYValue());
                            eV.add(xye.getError());
                        }
                    }
                    ExperimentData sdExpData = new ExperimentData(null, resonanceSource, xV, yV, eV);
                    experimentalDataSets.add(sdExpData);
                }
                var orderPars = atom.getOrderPars();
                parValues = new ArrayList<>();
                if (orderPars != null) {
                    for (var entry : orderPars.entrySet()) {
                        var orderPar = entry.getValue();
                        if (!chartInfo.mapName.contains(entry.getKey().name())) {
                            continue;
                        }
                        String modelName = orderPar.getModel();
                        var model = MFModelIso.buildModel(modelName, true, 0.0, 0.0, false);
                        var parNames = model.getParNames();
                        double[] pars = new double[parNames.size()];
                        double[] errs = new double[parNames.size()];
                        int iPar = 0;
                        for (var parName : parNames) {
                            pars[iPar] = orderPar.getValue(parName);
                            Double err = orderPar.getError(parName);
                            errs[iPar] = err == null ? 0.0 : err;
                            ParValue parValue = new ParValue(resonanceSource, "", parName, pars[iPar], errs[iPar]);
                            parValues.add(parValue);
                            iPar++;
                        }
                        Double aicc = orderPar.getAICC();
                        Double aic = orderPar.getAIC();
                        Double rms = orderPar.getRMS();
                        Double rChiSq = orderPar.getReducedChiSqr();
                        Integer nValues = orderPar.getN();
                        updateFitQuality(aic, aicc, rms, rChiSq, nValues);
                        double[] extras = new double[1];
                        var guiPlotEquation = new GUIPlotEquation(modelName, "spectralDensity", pars, errs, extras);
                        guiPlotEquation.setColor(PlotData.colors[iSeries % 8]);
                        equations.add(guiPlotEquation);
                        iSeries++;

                    }
                }
            }

        }
        chartInfo.setStates(allStates);
        updateTable(experimentalDataSets, fitMode);
        if (chartInfo.hasResidues()) {
            setControls();
            if (parValues != null) {
                xychart.setNames("Spectral Density", "\u03C9 (1/ns)", "log10[J(\u03C9)/1ns]", "0");
                updateTableWithPars(parValues);
            } else {
                updateTableWithPars(chartInfo);
            }
            updateEquation(chartInfo.mapName.get(0), chartInfo.getResidues(), chartInfo.equationName);
        }
        plotData.setData(allData);
        setBounds();
        plotData.autoScale(false);
        plotData.setEquations(equations);
    }

    public double getMaxY(ExperimentSet experimentSet, String equationName, String mapName, String state, ResonanceSource[] resSources) {
        double maxValue = Double.NEGATIVE_INFINITY;
        if ((experimentSet != null) && (resSources != null)) {
            for (Experiment expData : experimentSet.getExperimentData()) {
                if (!ExperimentSet.matchStateString(state, expData.getState())) {
                    continue;
                }
                String expName = expData.getName();
                for (ResonanceSource resNum : resSources) {

                    GUIPlotEquation equation = ChartUtil.getEquation(expData,
                            mapName, resNum, equationName, expData.getState(),
                            expData.getNucleusField());
                    if (equation != null) {
                        double minX = equation.getMinX();
                        double valueY = equation.calculate(minX);
                        if (valueY > maxValue) {
                            maxValue = valueY;
                        }
                    } else {
                        DataSeries series = ChartUtil.getMapData(mapName, expName, resNum);
                        double valueY = series.getValues().stream().mapToDouble(XYValue::getYValue).max().getAsDouble();
                        if (valueY > maxValue) {
                            maxValue = valueY;
                        }

                    }
                }
            }
        }
        return maxValue;
    }

    public void setSimData(EquationFitter equationFitter) {
        ArrayList<Double> xValues = getSimXData();
        ArrayList<Double> yValues = getSimYData();
        ArrayList<Double> errValues = getSimErrData();
        double[] extras = simControls.getExtras();
        int nX;
        ArrayList[] allXValues;
        ArrayList<Double> fieldValuesX = new ArrayList<>();
        if (equationFitter instanceof CESTFitter) {
            nX = 1;
            allXValues = new ArrayList[extras.length + nX + 1];
            allXValues[0] = xValues;
            allXValues[3] = fieldValuesX;
            for (int i = 0; i < yValues.size(); i++) {
                double fieldValue = CoMDPreferences.getRefField() * simControls.getNucleus().getFreqRatio();
                fieldValuesX.add(fieldValue);
            }
        } else {
            nX = 3;
            allXValues = new ArrayList[extras.length + nX];
            ArrayList<Double> fieldValuesH = new ArrayList<>();
            allXValues[0] = xValues;
            allXValues[1] = fieldValuesX;
            allXValues[2] = fieldValuesH;
            for (int i = 0; i < yValues.size(); i++) {
                double fieldValue = CoMDPreferences.getRefField() * simControls.getNucleus().getFreqRatio();
                fieldValuesX.add(fieldValue);
                fieldValuesH.add(CoMDPreferences.getRefField());
            }
        }
        for (int j = 0; j < extras.length; j++) {
            ArrayList<Double> xValuesEx = new ArrayList<>();
            for (int i = 0; i < yValues.size(); i++) {
                xValuesEx.add(extras[j]);
            }
            allXValues[nX + j] = xValuesEx;
        }
        equationFitter.setData(allXValues, yValues, errValues);
    }

    public void guessSimData() {
        EquationFitter equationFitter = getFitter();
        setSimData(equationFitter);
        String equationName = simControls.getEquation();
        List<ParValueInterface> guesses = equationFitter.guessPars(equationName);
        if (guesses != null) {
            simControls.updateSliders(guesses, equationName);
            simControls.simSliderAction("");
        }
    }

    public void fitSimData() {
        EquationFitter equationFitter = getFitter();
        setSimData(equationFitter);
        String equationName = simControls.getEquation();
        equationFitter.setupFit(equationName);
        int[][] map = equationFitter.getFitModel().getMap();
        double[] sliderGuesses = null;
        if (sliderGuessCheckBox.isSelected()) {
            sliderGuesses = simControls.sliderGuess(equationName, map);
        }
        CoMDOptions options = new CoMDOptions(true);
        fitResult = equationFitter.doFit(equationName, sliderGuesses, options);
        updateAfterFit(fitResult);
    }

    ArrayList<Double> getSimXData() {
        ObservableList<DataSeries> allData = xychart.getData();
        ArrayList<Double> xValues = new ArrayList<>();
        for (DataSeries series : allData) {
            for (XYValue dataPoint : series.getData()) {
                double x = dataPoint.getXValue();
                xValues.add(x);
            }
        }
        return xValues;
    }

    ArrayList<Double> getSimYData() {
        ObservableList<DataSeries> allData = xychart.getData();
        ArrayList<Double> yValues = new ArrayList<>();
        for (DataSeries series : allData) {
            for (XYValue dataPoint : series.getData()) {
                double y = dataPoint.getYValue();
                yValues.add(y);
            }
        }
        return yValues;
    }

    ArrayList<Double> getSimErrData() {
        ObservableList<DataSeries> allData = xychart.getData();
        ArrayList<Double> errValues = new ArrayList<>();
        for (DataSeries series : allData) {
            for (XYValue dataPoint : series.getData()) {
                double errValue = 0.0;
                if (dataPoint instanceof XYEValue) {
                    XYEValue xyeValue = (XYEValue) dataPoint;
                    errValue = xyeValue.getError();
                }
                errValues.add(errValue);
            }
        }
        return errValues;
    }

    @FXML
    void showSimData() {
        ObservableList<DataSeries> allData = FXCollections.observableArrayList();
        String equationName = simControls.getEquation();
        String simMode = getSimMode();
        EquationType eType = ResidueFitter.getEquationType(simMode, equationName);
        int[][] map = eType.makeMap(1);
        EquationFitter equationFitter = getFitter();
        equationFitter.getFitModel().setMap(map);
        double[] sliderGuesses = simControls.sliderGuess(equationName, map);
        double[] yBounds = xychart.getYBounds();
        double sdev = Math.abs(yBounds[1] - yBounds[0]) * 0.02;
        if (genDataSDevTextField.getText().equals("")) {
            genDataSDevTextField.setText(String.valueOf(sdev));
        }
        double[] xValues;
        if (simMode.equals("cest") || simMode.equals("r1rho")) {
            int nPts = Integer.parseInt(genDataNPtsTextField.getText());
            double xLB = Double.parseDouble(genDataXLBTextField.getText());
            double xUB = Double.parseDouble(genDataXUBTextField.getText());
            xValues = equationFitter.getSimX(nPts, xLB, xUB);
        } else {
            String[] xvals = genDataXValTextField.getText().split(" ");
            double[] xVals = new double[xvals.length];
            for (int i = 0; i < xvals.length; i++) {
                xVals[i] = Double.parseDouble(xvals[i]);
            }
            xValues = xVals;
        }
        double fieldRef;

        List<DataSeries> data = new ArrayList<>();
        DataSeries series = new DataSeries();
        series.setName("sim" + ":" + "0");
        data.add(series);
        int nX = 3;
        for (PlotEquation eqn : xychart.plotEquations) {
            double[] extras = eqn.getExtras();
            double[] ax = new double[1 + extras.length];
            if (extras.length - 1 >= 0) System.arraycopy(extras, 0, ax, 1, extras.length);
            for (double xValue : xValues) {
                ax[0] = xValue;
                double yValue = eqn.calculate(sliderGuesses, ax);
                yValue += Double.parseDouble(genDataSDevTextField.getText()) * rand.nextGaussian(); //sdev * rand.nextGaussian();
                XYValue dataPoint = new XYEValue(xValue, yValue, Double.parseDouble(genDataSDevTextField.getText()));
                dataPoint.setExtraValue(sdev);
                series.add(dataPoint);
            }
        }
        allData.addAll(data);
        xychart.setData(allData);
    }

    @FXML
    public void loadSimData() {
        ObservableList<DataSeries> allData = FXCollections.observableArrayList();
        FileChooser fileChooser = new FileChooser();
        File file = fileChooser.showOpenDialog(null);
        if (file != null) {
            try {
                List<Double>[] dataValues = DataIO.loadSimDataFile(file);
                List<DataSeries> data = new ArrayList<>();
                DataSeries series = new DataSeries();
                series.setName("sim" + ":" + "0");
                data.add(series);
                for (int i = 0; i < dataValues[0].size(); i++) {
                    XYValue dataPoint;
                    if (!dataValues[2].isEmpty()) {
                        dataPoint = new XYEValue(dataValues[0].get(i),
                                dataValues[1].get(i), dataValues[2].get(i));
                    } else {
                        dataPoint = new XYValue(dataValues[0].get(i), dataValues[1].get(i));

                    }
                    series.add(dataPoint);
                }
                allData.addAll(data);
                xychart.setData(allData);
            } catch (IOException ioE) {
                ExceptionDialog dialog = new ExceptionDialog(ioE);
                dialog.showAndWait();
            }
        }
    }

    @FXML
    public void showMCplot() {
        if (bootstrapSamplePlots == null) {
            bootstrapSamplePlots = new BootstrapSamplePlots(this);
        }
        bootstrapSamplePlots.showMCplot();
    }

    @FXML
    private void showPreferences() {
        if (preferencesController == null) {
            preferencesController = PreferencesController.create(primaryStage);
        }
        if (preferencesController != null) {
            preferencesController.getStage().show();
        } else {
            System.out.println("Coudn't make controller");
        }
    }

    @FXML
    private void showConsole() {
        ConsoleController.getConsoleController().show();
    }

    @FXML
    private void exitProgram() {
        MainApp.checkExit();
    }

    @FXML
    private void printXYChart() {
        printChart(chartPane);
    }

    @FXML
    private void printBarChart() {
        printChart(chartBox);
    }

    private void printChart(Node node) {
        PrinterJob job = PrinterJob.createPrinterJob();
        if (job != null) {
            boolean selected = job.showPrintDialog(null);
            if (selected) {
                boolean pageSetup = job.showPageSetupDialog(null);
                if (pageSetup) {
                    boolean success = job.printPage(node);
                    if (success) {
                        job.endJob();
                    }
                }
            }
        }
    }

    public File getInitialDirectory() {
        if (initialDir == null) {
            String homeDirName = System.getProperty("user.home");
            initialDir = new File(homeDirName);
        }
        return initialDir;
    }

    @FXML
    void exportSVGAction() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Export to SVG");
        fileChooser.setInitialDirectory(getInitialDirectory());
        File selectedFile = fileChooser.showSaveDialog(null);
        if (selectedFile != null) {
            SVGGraphicsContext svgGC = new SVGGraphicsContext();
            try {
                Canvas canvas = xychart.getCanvas();
                svgGC.create(canvas.getWidth(), canvas.getHeight(), selectedFile.toString());
                xychart.exportVectorGraphics(svgGC);
                svgGC.saveFile();
            } catch (GraphicsIOException ex) {
                ExceptionDialog eDialog = new ExceptionDialog(ex);
                eDialog.showAndWait();
            }
        }
    }

    @FXML
    void exportBarPlotSVGAction() {
        FileChooser fileChooser = new FileChooser();
        fileChooser.setTitle("Export to SVG");
        fileChooser.setInitialDirectory(getInitialDirectory());
        File selectedFile = fileChooser.showSaveDialog(null);
        if (selectedFile != null) {
            SVGGraphicsContext svgGC = new SVGGraphicsContext();
            try {
                svgGC.create(barPlotCanvas.getWidth(), barPlotCanvas.getHeight(), selectedFile.toString());
                exportBarPlotSVGAction(svgGC);
                svgGC.saveFile();
            } catch (GraphicsIOException ex) {
                ExceptionDialog eDialog = new ExceptionDialog(ex);
                eDialog.showAndWait();
            }
        }
    }

    protected void exportBarPlotSVGAction(SVGGraphicsContext svgGC) throws GraphicsIOException {
        svgGC.beginPath();
        for (ResidueChart resChart : barCharts) {
            resChart.drawChart(svgGC);
        }
        if (ssPainter != null) {
            double ssHeight = ssPainter.getHeight();
            Axis axis = activeChart.xAxis;
            ssPainter.paintSS(svgGC, axis.getXOrigin(), barPlotCanvas.getHeight() - ssHeight,
                    axis.getWidth(),
                    axis.getLowerBound(), axis.getUpperBound());
        }
    }

    public boolean hasExperimentSet() {
        return (chartInfo.valueSet != null) && (chartInfo.valueSet instanceof ExperimentSet);
    }

    /**
     * @return the valueSet
     */
    public ExperimentSet getCurrentExperimentSet() {
        if (chartInfo.valueSet instanceof ExperimentSet) {
            return (ExperimentSet) chartInfo.valueSet;
        } else {
            return null;
        }
    }

    /**
     * @param valueSet the valueSet to set
     */
    public void setCurrentExperimentSet(ValueSet valueSet) {
        this.chartInfo.valueSet = valueSet;
    }

    public void interpolateRates() {
        String fieldStr = GUIUtils.input("Target Field");
        if (!fieldStr.isBlank()) {
            double field = Double.valueOf(fieldStr);
            OrderParameterTool.interpolateRates(field);
        }
    }

    public void calculateSpectralDensities() {
        OrderParameterTool.calculateSpectralDensities();
    }

}
