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

import org.comdnmr.modelfree.RelaxEquations;
import org.nmrfx.utils.properties.DirectoryOperationItem;
import org.nmrfx.utils.properties.ChoiceOperationItem;
import org.nmrfx.utils.properties.IntRangeOperationItem;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.prefs.Preferences;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.scene.Scene;
import javafx.scene.layout.Pane;
import javafx.stage.Stage;
import javafx.stage.StageStyle;
import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.data.Experiment;
import org.controlsfx.control.PropertySheet;
import org.nmrfx.utils.properties.DoubleRangeOperationItem;
import org.nmrfx.utils.properties.NvFxPropertyEditorFactory;
import org.nmrfx.utils.properties.BooleanOperationItem;

/**
 *
 * @author johnsonb
 */
public class PreferencesController implements Initializable {

    @FXML
    PropertySheet prefSheet;
    ChangeListener<String> stringListener;
    ChangeListener<String> datasetListener;
    ChangeListener<String> optimizerListener;
    ChangeListener<String> bootStrapOptimizerListener;
    ChangeListener<Integer> nprocessListener;
    ChangeListener<Boolean> cestEqnListener;
    ChangeListener<Boolean> cpmgEqnListener;
    ChangeListener<Boolean> r1rhoEqnListener;
    ChangeListener<Boolean> expEqnListener;
    Stage stage;

    static File nestaNMR = null;
    static File datasetDir = null;
    private static Map<String, String> recentMap = new HashMap<>();
    static String location = null;
    static Integer nProcesses = null;
    DoubleRangeOperationItem maxFreqItem;

    @Override
    public void initialize(URL url, ResourceBundle rb) {
        prefSheet.setPropertyEditorFactory(new NvFxPropertyEditorFactory());
        prefSheet.setMode(PropertySheet.Mode.CATEGORY);
        prefSheet.setModeSwitcherVisible(false);
        prefSheet.setSearchBoxVisible(false);

        datasetListener = (ObservableValue<? extends String> observableValue, String string, String string2) -> {
            setDatasetDirectory(new File(string2.trim()));
        };
        optimizerListener = (ObservableValue<? extends String> observableValue, String string, String string2) -> {
            CoMDPreferences.setOptimizer(string2.trim());
        };
        bootStrapOptimizerListener = (ObservableValue<? extends String> observableValue, String string, String string2) -> {
            CoMDPreferences.setBootStrapOptimizer(string2.trim());
        };
        nprocessListener = (ObservableValue<? extends Integer> observableValue, Integer n1, Integer n2) -> {
            setNProcesses(n2);
        };
        cestEqnListener = (ObservableValue<? extends Boolean> observableValue, Boolean cest1, Boolean cest2) -> {
            BooleanOperationItem item = (BooleanOperationItem) observableValue;
            CoMDPreferences.setCESTEquationState(item.getName(), cest2);
            PyController.mainController.updateEquationChoices();
        };
        cpmgEqnListener = (ObservableValue<? extends Boolean> observableValue, Boolean cpmg1, Boolean cpmg2) -> {
            BooleanOperationItem item = (BooleanOperationItem) observableValue;
            CoMDPreferences.setCPMGEquationState(item.getName(), cpmg2);
            PyController.mainController.updateEquationChoices();
        };
        r1rhoEqnListener = (ObservableValue<? extends Boolean> observableValue, Boolean r1rho1, Boolean r1rho2) -> {
            BooleanOperationItem item = (BooleanOperationItem) observableValue;
            CoMDPreferences.setR1RhoEquationState(item.getName(), r1rho2);
            PyController.mainController.updateEquationChoices();
        };
        expEqnListener = (ObservableValue<? extends Boolean> observableValue, Boolean exp1, Boolean exp2) -> {
            BooleanOperationItem item = (BooleanOperationItem) observableValue;
            CoMDPreferences.setExpEquationState(item.getName(), exp2);
            PyController.mainController.updateEquationChoices();
        };
        ArrayList<String> optimizers = new ArrayList<>();
        optimizers.add("CMA-ES");
        optimizers.add("BOBYQA");
        maxFreqItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            CoMDPreferences.setCPMGMaxFreq((Double) newV);
        }, CoMDPreferences.getCPMGMaxFreq(), 100.0, 5000.0, 100.0, 20000.0, "Limits", "CPMG Max Freq", "Max Frequency");
        DoubleRangeOperationItem rexRatioItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            CoMDPreferences.setRexRatio((Double) newV);
        }, CoMDPreferences.getRexRatio(), 0.0, 10.0, "Limits", "CPMG Rex Ratio", "Rex must be this many times rmsd");
        DoubleRangeOperationItem deltaABdiffItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            CoMDPreferences.setDeltaABDiff((Double) newV);
        }, CoMDPreferences.getDeltaABDiff(), 0.0, 2.0, "Limits", "CEST/R1Rho DeltaAB Min Diff", "DeltaAB difference magnitude must be larger than this number for valid exchange");
        DoubleRangeOperationItem refFieldItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            CoMDPreferences.setRefField((Double) newV);
        }, CoMDPreferences.getRefField(), 100.0, 1200.0, "Limits", "Reference Field", "Rex value reported at this H1 field");
        IntRangeOperationItem nSamplesItem = new IntRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            CoMDPreferences.setSampleSize((Integer) newV);
        }, CoMDPreferences.getSampleSize(), 10, 500, "Optimizer", "Bootstrap Samples", "Number of bootstrap samples");
        ChoiceOperationItem optimizerChoiceItem = new ChoiceOperationItem(prefSheet, optimizerListener,
                CoMDPreferences.getOptimizer(), optimizers, "Optimizer", "Refine Optimizer", "Optimizer for refinement");
        ChoiceOperationItem bootStrapOptimizerrChoiceItem = new ChoiceOperationItem(prefSheet, bootStrapOptimizerListener,
                CoMDPreferences.getBootStrapOptimizer(), optimizers, "Optimizer", "Bootstrap Optimizer", "Optimizer for bootstrap");

        DirectoryOperationItem locationFileItem = new DirectoryOperationItem(prefSheet, datasetListener, getDatasetDirectory().getPath(), "File Locations", "Datasets", "desc");

        int nProcessesDefault = Runtime.getRuntime().availableProcessors() / 2;
        IntRangeOperationItem nProcessesItem = new IntRangeOperationItem(prefSheet, nprocessListener, nProcessesDefault, 1, 32, "Processor", "NProcesses", "How many parallel processes to run during processing");

        ArrayList<String> cestEqnChoices = new ArrayList<>();
        cestEqnChoices.addAll(Arrays.asList("NOEX", "TROTT_PALMER", "SD", "BALDWINKAY", "LAGUERRE",
                "EIGENEXACT1", "EXACT0", "EXACT1", "EXACT2"));

        ArrayList<String> cpmgEqnChoices = new ArrayList<>();
        cpmgEqnChoices.addAll(Arrays.asList("NOEX", "CPMGFAST", "CPMGSLOW", "CPMGMQ"));

        ArrayList<String> r1rhoEqnChoices = new ArrayList<>();
        r1rhoEqnChoices.addAll(Arrays.asList("NOEX", "TROTT_PALMER", "BALDWINKAY", "LAGUERRE", "EXACT", "EXACT0"));

        ArrayList<String> expEqnChoices = new ArrayList<>();
        expEqnChoices.addAll(Arrays.asList("EXPAB", "EXPABC"));
//        prefSheet.getItems().addAll(locationTypeItem, locationFileItem, nProcessesItem, maxFreqItem, rexRatioItem, nSamplesItem);

        DoubleRangeOperationItem startingRadiusItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            CoMDPreferences.setStartingRadius((Double) newV);
        }, CoMDPreferences.getStartRadius(), 1.0, 30.0, "Optimizer", "Starting Radius", "Iniital search radius for optimizer");
        DoubleRangeOperationItem finalRadiusItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            CoMDPreferences.setFinalRadius((Double) newV);
        }, CoMDPreferences.getFinalRadius(), -6, -1.0, "Optimizer", "Final Radius", "Final search radius for BOBYQA");
        DoubleRangeOperationItem toleranceItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            CoMDPreferences.setTolerance((Double) newV);
        }, CoMDPreferences.getTolerance(), -7.0, -1.0, "Optimizer", "Tolerance", "Final tolerance for CMA-ES");

        BooleanOperationItem weightFitItem = new BooleanOperationItem(prefSheet,
                (obs, oldV, newV) -> {
                    CoMDPreferences.setWeightFit((Boolean) newV);
                }, CoMDPreferences.getWeightFit(), "Optimizer", "Weight Fit",
                "Weight residuals by error values");
        BooleanOperationItem nonParametricItem = new BooleanOperationItem(prefSheet,
                (obs, oldV, newV) -> {
                    CoMDPreferences.setNonParametric((Boolean) newV);
                }, CoMDPreferences.getNonParametricBootstrap(), "Optimizer", "Non-Parametric",
                "Bootstrapping with non-parametric method");
        BooleanOperationItem absValueItem = new BooleanOperationItem(prefSheet,
                (obs, oldV, newV) -> {
                    CoMDPreferences.setAbsValueFit((Boolean) newV);
                }, CoMDPreferences.getAbsValueFit(), "Optimizer", "Absolute Value",
                "Fit absolute value of deviations (not squared)");
        BooleanOperationItem neuralNetworkGuessItem = new BooleanOperationItem(prefSheet,
                (obs, oldV, newV) -> {
                    CoMDPreferences.setNeuralNetworkGuess((Boolean) newV);
                }, CoMDPreferences.getNeuralNetworkGuess(), "Optimizer", "Neural Network Guess",
                "Use Neural Network for guesses");
        BooleanOperationItem corrR1RhoItem = new BooleanOperationItem(prefSheet,
                (obs, oldV, newV) -> {
                    CoMDPreferences.setCalR1rhoCorr((Boolean) newV);
                }, CoMDPreferences.getCalR1rhoCorr(), "R1Rho Equations", "Correct R1rho",
                "Correction if starting in ground state");

        DoubleRangeOperationItem hnDistanceItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            RelaxEquations.setR("H", "N", (Double) newV * 1.0e-10);
        }, RelaxEquations.getR("H", "N") * 1.0e10, 1.0, 1.04, "ModelFree", "HN Distance", "H-N bond distance");

        DoubleRangeOperationItem hcDistanceItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            RelaxEquations.setR("H", "C", (Double) newV * 1.0e-10);
        }, RelaxEquations.getR("H", "C") * 1.0e10, 1.04, 1.12, "ModelFree", "HC Distance", "H-N bond distance");

        DoubleRangeOperationItem nSigmaItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            RelaxEquations.setSigma( "N", (Double) newV * -1.0e-6);
        }, RelaxEquations.getSigma("N") * -1.0e6, 10, 220, "ModelFree", "N CSA", "N Chemical Shift Anisotropy (* -1.0)");

        DoubleRangeOperationItem cSigmaItem = new DoubleRangeOperationItem(prefSheet, (obs, oldV, newV) -> {
            RelaxEquations.setSigma("C", (Double) newV * -1.0e-6);
        }, RelaxEquations.getSigma("C") * -1.0e6, 10, 220, "ModelFree", "C CSA", "C Chemical Shift Anisotropy (* -1.0)");

        prefSheet.getItems().addAll(nProcessesItem, refFieldItem, maxFreqItem, rexRatioItem, deltaABdiffItem,
                absValueItem, nonParametricItem, nSamplesItem, optimizerChoiceItem, bootStrapOptimizerrChoiceItem,
                startingRadiusItem, toleranceItem, finalRadiusItem, weightFitItem, neuralNetworkGuessItem, corrR1RhoItem, hnDistanceItem, hcDistanceItem, nSigmaItem, cSigmaItem);
        for (String eqn : cestEqnChoices) {
            boolean defaultState = CoMDPreferences.getCESTEquationState(eqn);
            BooleanOperationItem cestEqnListItem = new BooleanOperationItem(prefSheet, cestEqnListener, defaultState, "CEST Equations", eqn, "List of equations to use during CEST Fitting");
            prefSheet.getItems().add(cestEqnListItem);
        }
        for (String eqn : cpmgEqnChoices) {
            boolean defaultState = CoMDPreferences.getCPMGEquationState(eqn);
            BooleanOperationItem cpmgEqnListItem = new BooleanOperationItem(prefSheet, cpmgEqnListener, defaultState, "CPMG Equations", eqn, "List of equations to use during CPMG Fitting");
            prefSheet.getItems().add(cpmgEqnListItem);
        }
        for (String eqn1 : r1rhoEqnChoices) {
            boolean defaultState = CoMDPreferences.getR1RhoEquationState(eqn1);
            BooleanOperationItem r1rhoEqnListItem = new BooleanOperationItem(prefSheet, r1rhoEqnListener, defaultState, "R1Rho Equations", eqn1, "List of equations to use during R1Rho Fitting");
            prefSheet.getItems().add(r1rhoEqnListItem);
        }

        for (String eqn1 : expEqnChoices) {
            boolean defaultState = CoMDPreferences.getExpEquationState(eqn1);
            BooleanOperationItem expEqnListItem = new BooleanOperationItem(prefSheet, expEqnListener, defaultState, "Exp Equations", eqn1, "List of equations to use during Exp Fitting");
            prefSheet.getItems().add(expEqnListItem);
        }

    }

    public static PreferencesController create(Stage parent) {
        FXMLLoader loader = new FXMLLoader(PyController.class.getResource("/fxml/RINGPreferencesScene.fxml"));
        PreferencesController controller = null;
        Stage stage = new Stage(StageStyle.DECORATED);

        try {
            Scene scene = new Scene((Pane) loader.load());
            stage.setScene(scene);
            scene.getStylesheets().add("/styles/Styles.css");

            controller = loader.<PreferencesController>getController();
            stage.setTitle("Preferences");

            stage.initOwner(parent);
            controller.stage = stage;
            stage.show();
        } catch (IOException ioE) {
            ioE.printStackTrace();
            System.out.println(ioE.getMessage());
        }

        return controller;
    }

    public Stage getStage() {
        return stage;
    }

    @FXML
    private void closeAction(ActionEvent event) {
        stage.close();
    }

    /**
     * Returns the Directory for datasets,
     *
     * @return
     */
    public static File getDatasetDirectory() {
        if (datasetDir == null) {
            Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
            String filePath = prefs.get("DATASET-DIR", null);
            if (filePath != null) {
                datasetDir = new File(filePath);
            } else {
                datasetDir = new File("");
            }
        }
        return datasetDir;
    }

    /**
     * Returns the Directory for datasets,
     *
     * @param file the file or null to remove the path
     */
    public static void setDatasetDirectory(File file) {
        Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
        if (file != null) {
            datasetDir = new File(file.getPath());
            prefs.put("DATASET-DIR", datasetDir.getPath());
        } else {
            datasetDir = null;
            prefs.remove("DATASET-DIR");
        }

    }

    /**
     * Returns the Directory for datasets,
     *
     * @return
     */
    public static String getLocation() {
        if (location == null) {
            Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
            String value = prefs.get("LOCATION-TYPE", null);
            if (value != null) {
                location = value;
            } else {
                location = "FID Directory";
            }
        }
        return location;
    }

    /**
     * Returns the Directory for datasets,
     *
     * @param value the file or null to remove the path
     */
    public static void setLocation(String value) {
        Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
        if (value != null) {
            location = value;
            prefs.put("LOCATION-TYPE", value);
        } else {
            location = null;
            prefs.remove("LOCATION-TYPE");
        }

    }

    /**
     * Saves recently opened datasets. The path is persisted in the OS specific
     * registry.
     *
     */
    public static List<Path> getRecentDatasets() {
        return getRecentFileItem("RECENT-DATASETS");
    }

    public static List<Path> getRecentProjects() {
        return getRecentFileItem("RECENT-PROJECTS");
    }

    public static void saveRecentProjects(String fileName) {
        saveRecentFileItems(fileName, "RECENT-PROJECTS");
    }

    public static void saveRecentDatasets(String fileName) {
        saveRecentFileItems(fileName, "RECENT-DATASETS");
    }

    public static void saveRecentFileItems(String fileName, String mode) {
        Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
        String recentFileString = recentMap.get(mode);
        if (recentFileString == null) {
            recentFileString = prefs.get(mode, "");
            recentMap.put(mode, recentFileString);
        }
        String[] recentDatasets = recentFileString.split("\n");
        Map<String, Long> datasetMap = new HashMap<>();
        for (String recentDatasetEntry : recentDatasets) {
            String[] entry = recentDatasetEntry.split(";");
            File file = new File(entry[0]);
            if (file.exists()) {
                datasetMap.put(entry[0], Long.valueOf(entry[1]));
            }
        }
        datasetMap.put(fileName, System.currentTimeMillis());
        StringBuilder sBuilder = new StringBuilder();
        datasetMap.entrySet().stream().sorted((e1, e2) -> Long.compare(e2.getValue(), e1.getValue())).limit(15).forEach(e1 -> {
            sBuilder.append(e1.getKey());
            sBuilder.append(';');
            sBuilder.append(String.valueOf(e1.getValue()));
            sBuilder.append("\n");
        });
        recentFileString = sBuilder.toString();
        recentMap.put(mode, recentFileString);
        prefs.put(mode, recentFileString);
    }

    public static List<Path> getRecentFileItem(String mode) {
        String recentFileString = recentMap.get(mode);
        if (recentFileString == null) {
            Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
            recentFileString = prefs.get(mode, "");
            recentMap.put(mode, recentFileString);
        }
        String[] recentDatasets = recentFileString.split("\n");
        List<Path> result = new ArrayList<>();
        for (String recentDatasetEntry : recentDatasets) {
            String[] entry = recentDatasetEntry.split(";");
            File file = new File(entry[0]);
            if (file.exists()) {
                Path path = file.toPath();
                result.add(path);
            }
        }
        return result;
    }

    /**
     * Returns the Directory for datasets,
     *
     * @return
     */
    public static Integer getNProcesses() {
        if (nProcesses == null) {
            Preferences prefs = Preferences.userNodeForPackage(MainApp.class);
            String value = prefs.get("NPROCESSES", null);
            if (value != null) {
                nProcesses = Integer.parseInt(value);
            } else {
                nProcesses = Runtime.getRuntime().availableProcessors() / 2;
            }
        }
        return nProcesses;
    }

    /**
     * Returns the Directory for datasets,
     *
     * @param value the number of processes to use
     */
    public static void setNProcesses(Integer value) {
        Preferences prefs = Preferences.userNodeForPackage(Experiment.class);
        if (value != null) {
            nProcesses = value;
            prefs.put("NPROCESSES", String.valueOf(value));
        } else {
            nProcesses = null;
            prefs.remove("NPROCESSES");
        }

    }

}
