# -*- coding: utf-8 -*-
"""
Helper functions for script to generate Figure 3 of:

    Lithospheric anisotropy characterises the post-subduction setting of
    northern Borneo: new results from XKS splitting analysis
    Bacon, C.A., Rawlinson, N., Pilia, S., Gilligan, A., Wehner, D.,
    Cornwell, D.G., Tongkul, F.

which has been submitted to Earth and Planetary Science Letters.

"""

import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd


def plot_results_scatter(ax, df, clr, category, shape):
    """
    Plot the splitting intensity measurements.
    """
    
    # Get error bars
    pos_err = df["SplitInt_err_max"].values - df["SplitInt"].values
    neg_err = df["SplitInt"].values - df["SplitInt_err_min"].values
    yerrs = np.array([neg_err, pos_err])

    ax.errorbar(df["baz"].values, df["SplitInt"].values, yerr=yerrs,
                fmt=shape, capsize=5, color=clr, linewidth=1.5, ms=5,
                label=f"{category.title()} results")


def plot_phi_baz(ax, df, clr, category, shape):
    """
    Plot the splitting intensity measurements.
    """
    phis, bazs, pos_err, neg_err = [], [], [], []
    
    for i, row in df.iterrows():
        phis.append(row["phi"])
        bazs.append(row["baz"])
        pos_err.append(row["err_phi_max"] - row["phi"])
        neg_err.append(row["phi"] - row["err_phi_min"])
        if row["err_phi_max"] > 180:
            excess = row["err_phi_max"] - 180
            phis.append(-5)
            bazs.append(row["baz"])
            pos_err.append(excess+5)
            neg_err.append(0)
        if row["err_phi_min"] < 0:
            excess = row["err_phi_min"] + 0
            phis.append(185)
            bazs.append(row["baz"])
            pos_err.append(0)
            neg_err.append(abs(excess-5))
    yerrs = np.array([neg_err, pos_err])

    ax.errorbar(bazs, phis, yerr=yerrs,
                fmt=shape, capsize=5, ms=5, color=clr, linewidth=1.5,
                label=f"{category.title()} results")
