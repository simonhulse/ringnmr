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
package org.comdnmr.data;

import java.util.ArrayList;
import java.util.List;
import org.comdnmr.util.DataUtil;
import org.nmrfx.chemistry.Entity;
import org.nmrfx.chemistry.Residue;
import org.nmrfx.chemistry.relax.ResonanceSource;

/**
 *
 * @author Bruce Johnson
 */
public class ExperimentData {

    Experiment experiment;
    ResonanceSource dynSource;
    double[][] xValues;
    double[] errValues;
    double[] yValues;

    public ExperimentData(Experiment experiment, ResonanceSource dynSource,
            double[][] x, double[] y, double[] err) {
        this.experiment = experiment;
        this.dynSource = dynSource;
        this.xValues = DataUtil.clone2DArray(x);
        this.yValues = y.clone();
        this.errValues = err.clone();
    }

    public ExperimentData(Experiment experiment, ResonanceSource dynSource,List<DataIO.XYErrValue> xyErrValueList) {
        this.experiment = experiment;
        this.dynSource = dynSource;
        int nValues = xyErrValueList.size();
        this.xValues = new double[1][nValues];
        this.yValues = new double[nValues];
        this.errValues = new double[nValues];
        for (int i = 0; i < nValues; i++) {
            xValues[0][i] = xyErrValueList.get(i).x();
            yValues[i] = xyErrValueList.get(i).y();
            errValues[i] = xyErrValueList.get(i).err();
        }
    }

    public ExperimentData(Experiment experiment, ResonanceSource dynSource, List<DataIO.XArrayYErrValue> xArrayYErrValues, int nX) {
        this.experiment = experiment;
        this.dynSource = dynSource;
        int nValues = xArrayYErrValues.size();
        this.xValues = new double[nX][nValues];
        this.yValues = new double[nValues];
        this.errValues = new double[nValues];
        for (int i = 0; i < nValues; i++) {
            DataIO.XArrayYErrValue xArrayYErrValue = xArrayYErrValues.get(i);
            for (int j = 0; j < nX; j++) {
                xValues[j][i] = xArrayYErrValue.x()[j];
            }
            yValues[i] = xArrayYErrValue.y();
            errValues[i] = xArrayYErrValue.err();
        }
    }

    public ResonanceSource getSource() {
        return dynSource;
    }

    public double[][] getXValues() {
        return xValues;
    }

    public double[] getYValues() {
        return yValues;
    }

    public double[] getErrValues() {
        return errValues;
    }
    
    public void setErrValue(int i, double value) {
        errValues[i] = value;
    }

    public class DataValue {

        int index;
        ExperimentData resInfo;

        public DataValue(ExperimentData resInfo, int index) {
            this.resInfo = resInfo;
            this.index = index;
        }

        public int getIndex() {
            return index;
        }

        public double getX0() {
            return resInfo.xValues[0][index];
        }

        public double getX1() {
            double x1 = resInfo.xValues[0][index];
            if (resInfo.xValues.length > 1) {
                x1 = resInfo.xValues[1][index];
            }
            return x1;
        }

        public double getX2() {
            double x2 = resInfo.xValues[0][index];
            if (resInfo.xValues.length > 2) {
                x2 = resInfo.xValues[2][index];
            }
            return x2;
        }

        public double getY() {
            return resInfo.yValues[index];
        }

        public double getError() {
            return resInfo.errValues[index];
        }

        public String getRefPeak() {
            return resInfo.dynSource.getPeak() == null ? "" : resInfo.dynSource.getPeak().getName();
        }

        public int getPeak() {
            return resInfo.dynSource.getPeak() == null ? -1 : resInfo.dynSource.getPeak().getIdNum();
        }

        public String getName() {
            return resInfo.experiment == null ? "" : resInfo.getDataName();
        }

        public String getResidue() {
            if (resInfo == null) {
                return "1";
            } else {
                Entity entity = resInfo.dynSource.getAtoms()[0].getEntity();
                if (entity instanceof Residue) {
                    Residue residue = (Residue) entity;
                    return residue.getNumber();
                } else {
                    return "1";
                }
            }
        }

        public String getResName() {
            return resInfo == null ? "" : resInfo.dynSource.getAtoms()[0].getResidueName();
        }

        public String getAtomName() {
            return resInfo == null ? "" : resInfo.dynSource.getAtoms()[0].getName();
        }

        public ResonanceSource getResonanceSource() {
            return resInfo == null ? null : resInfo.getSource();
        }
    }

    public ArrayList<DataValue> getDataValues() {
        ArrayList<DataValue> dataValues = new ArrayList<>();
        if (xValues != null) {
            for (int i = 0; i < xValues[0].length; i++) {
                dataValues.add(new DataValue(this, i));
            }
        }
        return dataValues;
    }

    public String getDataName() {
        return experiment.getName();
    }

    public Experiment getExperimentData() {
        return experiment;
    }

    public ResonanceSource getDynSource() {
        return dynSource;
    }

}
