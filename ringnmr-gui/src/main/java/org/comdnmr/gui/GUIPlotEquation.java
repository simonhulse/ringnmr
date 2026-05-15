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

import javafx.scene.paint.Color;
import org.comdnmr.eqnfit.PlotEquation;

/**
 *
 * @author brucejohnson
 */
public class GUIPlotEquation extends PlotEquation {

    private Color color = Color.BLACK;
    private double scaleValue = 1.0;
    private double lineWidth = 1.0;

    public GUIPlotEquation(String expType, String name, double[] pars, double[] errs, double[] extras) {
        super(expType, name, pars, errs, extras);
    }

    public GUIPlotEquation(String expType, PlotEquation eqn) {
        super(expType, eqn.getName(), eqn.getPars(), eqn.getErrs(), eqn.getExtras());
    }

    /**
     * @return the color
     */
    public Color getColor() {
        return color;
    }

    /**
     * @param color the color to set
     */
    public void setColor(Color color) {
        this.color = color;
    }

    /**
     * @return the scaleValue
     */
    public double getScaleValue() {
        return scaleValue;
    }

    /**
     * @param scaleValue the scaleValue to set
     */
    public void setScaleValue(double scaleValue) {
        this.scaleValue = scaleValue;
    }

    public double getLineWidth() {
        return lineWidth;
    }

    public void setLineWidth(double lineWidth) {
        this.lineWidth = lineWidth;
    }

}
