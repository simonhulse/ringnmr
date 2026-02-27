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
package org.comdnmr.fit;

import org.comdnmr.data.DataIO;
import org.comdnmr.eqnfit.R1RhoFitter;
import org.comdnmr.eqnfit.CESTEquation;
import org.comdnmr.eqnfit.ExpEquation;
import org.comdnmr.eqnfit.EquationType;
import org.comdnmr.eqnfit.FitResult;
import org.comdnmr.eqnfit.R1RhoEquation;
import org.comdnmr.eqnfit.SSR1RhoEquation;
import org.comdnmr.eqnfit.CESTFitter;
import org.comdnmr.eqnfit.CurveFit;
import org.comdnmr.eqnfit.ExpFitter;
import org.comdnmr.eqnfit.EquationFitter;
import org.comdnmr.eqnfit.CPMGFitter;
import org.comdnmr.eqnfit.CPMGEquation;
import org.comdnmr.util.ProcessingStatus;
import org.comdnmr.data.ExperimentSet;
import org.comdnmr.data.ExperimentResult;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import javafx.beans.property.ReadOnlyObjectProperty;
import javafx.concurrent.Service;
import javafx.concurrent.Task;
import javafx.concurrent.Worker;
import org.comdnmr.data.Experiment;
import org.comdnmr.eqnfit.NOEEquation;
import org.comdnmr.util.CoMDOptions;
import org.nmrfx.chemistry.relax.RelaxTypes;
import org.nmrfx.chemistry.relax.ResonanceSource;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueFitter {

    final CoMDOptions options;
    private final FitResidues processDataset = new FitResidues();
    final ReadOnlyObjectProperty<Worker.State> stateProperty = processDataset.worker.stateProperty();
    private boolean isProcessing = false;
    ExperimentSet experimentSet;
    Function<Double, Double> updaterFunction;
    Function<ProcessingStatus, Double> statusFunction;
    List<List<ResonanceSource>> atomFitGroups = null;
    FitResult fitResult;

    public ResidueFitter(CoMDOptions options) {
        this.options = options;
    }

    public ResidueFitter(CoMDOptions options, Function<Double, Double> updaterFunction, Function<ProcessingStatus, Double> statusFunction) {
        this.updaterFunction = updaterFunction;
        this.statusFunction = statusFunction;
        this.options = options;
    }

    public void fitResiduesNow(ExperimentSet experimentSet) {
        this.experimentSet = experimentSet;
        experimentSet.setupMaps();
        this.atomFitGroups = null;
        experimentSet.clearResidueMap();
        fitAllAtoms(null);
        RelaxTypes relaxType = RelaxTypes.valueOf(experimentSet.getExpMode().toUpperCase());
        DataIO.addRelaxationFitResults(experimentSet, relaxType);
    }

    public void fitResidues(ExperimentSet experimentSet) {
        fitResidues(experimentSet, null);
    }

    public void fitResidues(ExperimentSet experimentSet, List<List<ResonanceSource>> residueFitGroups) {
        this.experimentSet = experimentSet;
        experimentSet.setupMaps();
        this.atomFitGroups = residueFitGroups;
        setProcessingOn();
        updateProgress(0.0);
        if (residueFitGroups == null) {
            experimentSet.clearResidueMap();
        }
        ((Service) processDataset.worker).restart();
    }

    public void updateProgress(Double f) {
        if (updaterFunction != null) {
            updaterFunction.apply(f);
        }
    }

    public void updateStatus(String s) {
        setProcessingStatus(s, true);
    }

    public void setProcessingStatus(String s, boolean ok) {
        setProcessingStatus(s, ok, null);
    }

    public void setProcessingStatus(String s, boolean ok, Throwable throwable) {
        if (statusFunction != null) {
            ProcessingStatus status = new ProcessingStatus(s, ok, throwable);
            statusFunction.apply(status);
        }
    }

    public void clearProcessingTextLabel() {
        ProcessingStatus status = new ProcessingStatus("");
        statusFunction.apply(status);
    }

    public void haltFit() {
        processDataset.worker.cancel();
    }

    synchronized void setProcessingOn() {
        isProcessing = true;
    }

    synchronized void setProcessingOff() {
        isProcessing = false;
    }

    boolean isProcessing() {
        return isProcessing;
    }

    void finishProcessing() {
        updateStatus("Done");
        try {
            String expMode = experimentSet.getExpMode().toUpperCase();
            if (expMode.equals("R1") || expMode.equals("R2")) {
                RelaxTypes relaxType = RelaxTypes.valueOf(experimentSet.getExpMode().toUpperCase());
                DataIO.addRelaxationFitResults(experimentSet, relaxType);
            }
        } catch (IllegalArgumentException iAE) {
            iAE.printStackTrace();

        }
    }

    public void fitAllResidueGroups(Task task) {
        for (List<ResonanceSource> dynFitGroup : atomFitGroups) {
            int nResidues = dynFitGroup.size();
            int nFit = 0;
            if (task.isCancelled()) {
                break;
            }
            ResonanceSource[] dynGroup = new ResonanceSource[dynFitGroup.size()];
            dynFitGroup.toArray(dynGroup);
            List<ExperimentResult> resInfoList = fitResidues(experimentSet, dynGroup, nFit, null);
            resInfoList.forEach((resInfo) -> {
                ResonanceSource dynSource = resInfo.getResonanceSource();
                if (!dynSource.deleted()) {
                    experimentSet.addExperimentResult(dynSource, resInfo);
                }
            });
            nFit++;
            updateProgress((1.0 * nFit) / nResidues);
        }
    }

    public List<List<ResonanceSource>> getAllAtoms() {
        Map<String, Experiment> expDataSets = experimentSet.getExperimentMap();
        Set<ResonanceSource> resSources = new TreeSet<>();
        expDataSets.values().forEach((expData) -> {
            expData.getDynamicsSources().stream().filter(dynSource -> !dynSource.deleted()).forEach(dynSource -> {
                resSources.add(dynSource);
            });
        });
        List<List<ResonanceSource>> allAtoms = new ArrayList<>();
        for (var resSource : resSources) {
            int groupIndex = -1;
            if (atomFitGroups != null) {
                int i = 0;
                for (List<ResonanceSource> atomGroup : atomFitGroups) {
                    if (atomGroup.contains(resSource)) {
                        groupIndex = i;
                        break;
                    }
                    i++;
                }
            }
            if (groupIndex != -1) {
                List<ResonanceSource> atomGroup = atomFitGroups.get(groupIndex);
                if (atomGroup.get(0).equals(resSource)) {
                    allAtoms.add(atomGroup);
                }
            } else if (atomFitGroups == null) {
                List<ResonanceSource> singleAtom = new ArrayList<>();
                singleAtom.add(resSource);
                allAtoms.add(singleAtom);
            }
        }

        return allAtoms;
    }

    public void fitAllAtoms(Task task) {
        List<List<ResonanceSource>> allAtoms = getAllAtoms();
        int nGroups = allAtoms.size();
        int nFit = 0;
        for (List<ResonanceSource> atomList : allAtoms) {
            if ((task != null) && task.isCancelled()) {
                break;
            }
            ResonanceSource[] atomGroup = new ResonanceSource[atomList.size()];
            atomList.toArray(atomGroup);
            List<ExperimentResult> resInfoList = fitResidues(experimentSet, atomGroup, nFit, null);
            resInfoList.forEach((resInfo) -> {
                ResonanceSource fitAtom = resInfo.getResonanceSource();
                experimentSet.addExperimentResult(fitAtom, resInfo);
            });
            nFit++;
            if (task != null) {
                updateProgress((1.0 * nFit) / nGroups);
            }
        }

    }

    EquationFitter getFitter(CoMDOptions options) {
        EquationFitter fitter;
        switch (experimentSet.getExpMode()) {
            case "cpmg":
                fitter = new CPMGFitter(options);
                break;
            case "r1":
            case "r2":
            case "rq":
            case "rap":
                fitter = new ExpFitter(options);
                break;
            case "cest":
                fitter = new CESTFitter(options);
                break;
            case "r1rho":
                fitter = new R1RhoFitter(options);
                break;
            default:
                throw new IllegalArgumentException("Invalid mode " + experimentSet.getExpMode());
        }
        return fitter;
    }

    public List<ExperimentResult> doNOE(ExperimentSet experimentSet, ResonanceSource[] dynSources, int groupId, String useEquation) {
        List<ExperimentResult> resInfoList = new ArrayList<>();
        Map<ResonanceSource, ExperimentResult> resMap = new HashMap<>();
        for (var dynSource : dynSources) {
            String expName = experimentSet.name();
            ExperimentResult residueInfo = new ExperimentResult(experimentSet, dynSource, groupId, dynSources.length);
            resInfoList.add(residueInfo);
            resMap.put(dynSource, residueInfo);
            residueInfo.addFitResult(null);
        }
        return resInfoList;

    }

    public List<ExperimentResult> fitResidues(ExperimentSet experimentSet, ResonanceSource[] dynSources, int groupId, String useEquation) {
        this.experimentSet = experimentSet;
        experimentSet.setupMaps();
        Map<String, FitResult> fitResults = new HashMap<>();
        double aicMin = Double.MAX_VALUE;
        String bestEquation = "NOEX";
        List<String> equationNames;
        switch (experimentSet.getExpMode()) {
            case "cpmg":
                equationNames = CPMGFitter.getEquationNames();
                break;
            case "r1":
            case "r2":
            case "rq":
            case "rap":
                equationNames = ExpFitter.getEquationNames();
                bestEquation = "EXPAB";
                break;
            case "cest":
                equationNames = CESTFitter.getEquationNames();
                bestEquation = "NOEX";
                break;
            case "r1rho":
                equationNames = R1RhoFitter.getEquationNames();
                bestEquation = "NOEX";
                break;
            case "noe":
                return doNOE(experimentSet, dynSources, groupId, useEquation);
            default:
                throw new IllegalArgumentException("Invalid mode " + experimentSet.getExpMode());
        }
        for (String equationName : equationNames) {
            if ((useEquation != null) && !equationName.equals(useEquation)) {
                continue;
            }

            EquationFitter equationFitter = getFitter(options);
            equationFitter.setData(experimentSet, dynSources);
            CoMDOptions options = new CoMDOptions(true);
            var fitResultOpt = equationFitter.doFit(equationName, null, options);
            if (fitResultOpt.isPresent()) {
                fitResult = fitResultOpt.get();
                fitResults.put(equationName, fitResult);
                if (fitResult.getAicc() < aicMin) {
                    aicMin = fitResult.getAicc();
                    if (fitResult.exchangeValid()) {
                        bestEquation = equationName;
                    }
                }
            }
        }
        List<ExperimentResult> resInfoList = new ArrayList<>();
        Map<ResonanceSource, ExperimentResult> resMap = new HashMap<>();
        for (var dynSource : dynSources) {
            ExperimentResult residueInfo = new ExperimentResult(experimentSet, dynSource, groupId, dynSources.length);
            resInfoList.add(residueInfo);
            resMap.put(dynSource, residueInfo);
            equationNames.forEach((equationName) -> {
                residueInfo.addFitResult(fitResults.get(equationName));
            });
        }

        for (String equationName : equationNames) {
            if ((useEquation != null) && !equationName.equals(useEquation)) {
                continue;
            }
            fitResult = fitResults.get(equationName);
            if (fitResult != null) {
                int nCurves = fitResult.getNCurves();
                for (int iCurve = 0; iCurve < nCurves; iCurve++) {
                    CurveFit curveFit = fitResult.getCurveFit(iCurve);
                    ResonanceSource dynSource = curveFit.getDynamicsSource();
                    ExperimentResult experimentResult = resMap.get(dynSource);
                    experimentResult.addCurveFit(curveFit, bestEquation.equals(equationName));
                }
            }
        }
        return resInfoList;
    }

    public FitResult getFitResult() {
        return fitResult;
    }

    private class FitResidues {

        String script;
        public Worker<Integer> worker;

        private FitResidues() {
            worker = new Service<Integer>() {

                @Override
                protected Task createTask() {
                    return new Task() {
                        @Override
                        protected Object call() {
                            updateStatus("Start processing");
                            updateTitle("Start Processing");
                            fitAllAtoms(this);
                            return 0;
                        }
                    };
                }
            };

            ((Service<Integer>) worker).setOnSucceeded(event -> {
                finishProcessing();
                setProcessingOff();
            });
            ((Service<Integer>) worker).setOnCancelled(event -> {
                setProcessingOff();
                setProcessingStatus("cancelled", false);
            });
            ((Service<Integer>) worker).setOnFailed(event -> {
                setProcessingOff();
                final Throwable exception = worker.getException();
                setProcessingStatus(exception.getMessage(), false, exception);

            });

        }
    }

    public static EquationType getEquationType(String expType, String name) throws IllegalArgumentException {
        EquationType equationType = null;
//        try {
        expType = expType.toLowerCase();
        switch (expType) {
            case "cpmg":
                equationType = CPMGEquation.valueOf(name);
                break;
            case "r1":
            case "t1":
                equationType = ExpEquation.valueOf(name);
                break;
            case "r2":
            case "t2":
                equationType = ExpEquation.valueOf(name);
                break;
            case "exp":
                equationType = ExpEquation.valueOf(name);
                break;
            case "cest":
                equationType = CESTEquation.valueOf(name);
                break;
            case "r1rho":
                equationType = R1RhoEquation.valueOf(name);
                break;
            case "ssr1rho":
                equationType = SSR1RhoEquation.valueOf(name);
                break;
            case "noe":
                equationType = NOEEquation.valueOf(name);
                break;
            default:
                break;
        }
//        } catch (IllegalArgumentException iaE) {
////            try {
////                equationType = ExpEquation.valueOf(name);
////            } catch (IllegalArgumentException iaE2) {
////                try {
////                    equationType = CESTEquation.valueOf(name);
////                } catch (IllegalArgumentException iaE3) {
////                    equationType = R1RhoEquation.valueOf(name);
////                }
////            }
//        }
        return equationType;
    }

}
