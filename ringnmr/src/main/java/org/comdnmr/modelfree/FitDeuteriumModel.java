package org.comdnmr.modelfree;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.nmrfx.chemistry.*;
import org.nmrfx.chemistry.relax.*;

import java.util.*;

public class FitDeuteriumModel extends FitModel {

    @Override
    protected DeuteriumStructureValues loadData() { return getData(false); }

    @Override
    protected String getNoDataMessage() {
        return "No relaxation data to analyze. Need R1, R2, RAP, and optionally RQ";
    }

    public static DeuteriumStructureValues getData(boolean requireCoords) {
        DeuteriumStructureValues molDataValues = new DeuteriumStructureValues();
        MoleculeBase mol = MoleculeFactory.getActive();
        if (mol != null) {
            for (Polymer polymer : mol.getPolymers()) {
                for (Residue residue : polymer.getResidues()) {
                    for (Atom atom : residue.getAtoms()) {
                        MolDataValues<? extends RelaxDataValue> molData = getMolDataValues(atom, requireCoords);
                        if ((molData != null) && !molData.getData().isEmpty()) {
                            molDataValues.put(molData.getSpecifier(), molData);
                        }
                    }
                }
            }
        }
        return molDataValues;
    }

    public static MolDataValues<? extends RelaxDataValue> getMolDataValues(Atom atom, boolean requireCoords) {
        var relaxData = atom.getRelaxationData();
        boolean hasDeuterium = relaxData.values().stream().
                anyMatch(v -> v.getRelaxationSet().relaxType() == RelaxTypes.RQ);
        MolDataValues<DeuteriumDataValue> molData = null;
        if (hasDeuterium) {
            Vector3D vec = null;
            if (requireCoords) {
                if (atom.getParent() != null) {
                    var aPoint = atom.getPoint();
                    var bPoint = atom.getParent().getPoint();
                    if ((aPoint == null) || (bPoint == null)) {
                        return null;
                    } else {
                        vec = aPoint.subtract(bPoint);
                    }
                } else {
                    return null;
                }
            }

            molData = vec == null
                    ? new DeuteriumMolDataValues(atom)
                    : new DeuteriumMolDataValues(atom, vec.toArray());
            Set<Long> fields = new HashSet<>();
            for (var entry : relaxData.entrySet()) {
                RelaxationSet relaxationSet = entry.getKey();
                if (relaxationSet.active()) {
                    fields.add(Math.round(relaxationSet.field()));
                }
            }
            for (var field : fields) {
                Double r1 = null;
                Double r1Error = null;

                Double r2 = null;
                Double r2Error = null;

                Double rQ = null;
                Double rQError = null;

                Double rAP = null;
                Double rAPError = null;

                for (var entry : relaxData.entrySet()) {
                    RelaxationData data = entry.getValue();
                    RelaxationSet relaxationSet = entry.getKey();
                    if (!relaxationSet.active()) {
                        continue;
                    }
                    if (Math.round(relaxationSet.field()) == field) {
                        switch (relaxationSet.relaxType()) {
                            case R1 -> {
                                r1 = data.getValue();
                                r1Error = data.getError();
                            }
                            case R2 -> {
                                r2 = data.getValue();
                                r2Error = data.getError();
                            }
                            case RQ -> {
                                rQ = data.getValue();
                                rQError = data.getError();
                            }
                            case RAP -> {
                                rAP = data.getValue();
                                rAPError = data.getError();
                            }
                        }
                    }
                }
                if ((r1 != null) && (r2 != null) && (rQ != null) && (rAP != null)) {
                    RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, "D", "C");

                    DeuteriumDataValue dValue = new DeuteriumDataValue(molData, r1, r1Error, r2, r2Error,
                            rQ, rQError, rAP, rAPError, relaxObj);
                    ((DeuteriumMolDataValues) molData).addData(dValue);
                }
            }
        }
        return molData;
    }
}
