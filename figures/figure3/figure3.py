# -*- coding: utf-8 -*-
"""
This script will create a summary figure for an example splitting measurement
Specifically, this script can be used to generate supplemental Figure S5 of:

    Lithospheric anisotropy characterises the post-subduction setting of
    northern Borneo: new results from XKS splitting analysis
    Bacon, C.A., Rawlinson, N., Pilia, S., Gilligan, A., Wehner, D.,
    Cornwell, D.G., Tongkul, F.

which has been submitted to Earth and Planetary Science Letters.

"""

import warnings
import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd

import functions as fig3_utils

try:
    plt.style.use("splitracer_summary")
except OSError:
    warnings.warn("You have not added the 'splitracer_summary.mplstyle' stylesheet to your"
                  " matplotlib config - continuing with matplotlib defaults.")
mpl.rcParams["font.family"] = "Helvetica"

# Read in data
results = pd.read_csv("files/results.txt", sep="\t", header=None)
energy_grid = pd.read_csv("files/energy_grid.txt", delim_whitespace=True, header=None)
dts, phis = energy_grid[0].values, energy_grid[1].values
Z = energy_grid[2].values.reshape(len(set(dts)), len(set(phis)))
diff_dt = (dts[-1] - dts[0]) / (len(dts) * 2)
diff_phi = (phis[-1] - phis[0]) / (len(phis) * 2)
X, Y = np.mgrid[dts[0]-diff_dt:dts[-1]+diff_dt:Z.shape[0]*1j,
                phis[0]-diff_phi:phis[-1]+diff_phi:Z.shape[1]*1j] 

phi, dt, minerr_phi, maxerr_phi, minerr_dt, maxerr_dt, _, _, conf_level = results[1].values

fast_err = np.mean([np.abs(phi - err) for err in [minerr_phi, maxerr_phi]])
delay_err = np.mean([np.abs(dt - err) for err in [minerr_dt, maxerr_dt]])

# Read in SI best-fitting results file
si_results = pd.read_csv("./files/SI_results.txt", header=None, delim_whitespace=True)
splitting_results = pd.read_csv("./files/splitting_results.txt", delim_whitespace=True)

# Split results by grades
good_results = splitting_results[splitting_results["category"] == "good"]
average_results = splitting_results[splitting_results["category"] == "average"]
poor_results = splitting_results[splitting_results["category"] == "poor"]
null_results = splitting_results[splitting_results["category"] == "null-measurement"]

splitting_results = splitting_results[splitting_results["category"] != "poor"]
splitting_results = splitting_results[splitting_results["err_phi_min"] != np.inf]
splitting_results = splitting_results[splitting_results["err_phi_max"] != np.inf]

splitting_results["phi"] = [phi - 180 if phi > 90 else phi
                            for phi in splitting_results["phi"].values]
splitting_results["err_phi_max"] = [phi + 180 if phi < 0 else phi
                                    for phi
                                    in splitting_results["err_phi_max"].values]
splitting_results["err_phi_min"] = [phi - 180 if phi > 90 else phi
                                    for phi
                                    in splitting_results["err_phi_min"].values]

good = splitting_results[splitting_results["category"] == "good"]
fair = splitting_results[splitting_results["category"] == "average"]
null = splitting_results[splitting_results["category"] == "null-measurement"]

fig = plt.figure(figsize=(16, 8), constrained_layout=True)
ax_dict = fig.subplot_mosaic(
    """
    AABB
    AACC
    """
)

ax = ax_dict["A"]
ax.pcolormesh(X, Y, Z, edgecolors="face", cmap="inferno_r", shading="gouraud")

ax.contour(X, Y, Z, 9, colors="#cccccc")
ax.contour(X, Y, Z, [conf_level,], colors="k")
ax.scatter(dt, phi, marker="+", s=200, label="Best-fitting (\u03C6, \u03B4t)")

ax.set_ylabel("\u03C6, °")
ax.set_xlabel("Delay time, s")
plt.xlim([0, 4])
plt.ylim([-90, 90])
ytickrange = np.arange(-90, 91, 15)
ax.set_yticks(ytickrange)
ax.set_yticklabels(ytickrange)

ax = ax_dict["C"]
# --- Plot splitting intensity best fit ---
# Unpack variables from file
fast_direction, delay_time, fast_err, delay_err = si_results.loc[0]

# Create splitting intensity function and uncertainties
x = np.arange(0, 360, 0.01)
best_fit = delay_time * np.sin(np.pi * (x - fast_direction) / 90)
fit_min1 = (delay_time - delay_err) * np.sin(np.pi * (x - (fast_direction - fast_err)) / 90)
fit_max1 = (delay_time + delay_err) * np.sin(np.pi * (x - (fast_direction + fast_err)) / 90)
fit_min2 = (delay_time - delay_err) * np.sin(np.pi * (x - (fast_direction + fast_err)) / 90)
fit_max2 = (delay_time + delay_err) * np.sin(np.pi * (x - (fast_direction - fast_err)) / 90)

fit_min, fit_max = [], []
for vals in zip(fit_min1, fit_min2, fit_max1, fit_max2):
    fit_min.append(min(vals))
    fit_max.append(max(vals))

# Plot
ax.plot(x, best_fit, label="SVD", color="black", linewidth=1.5, linestyle="--")
ax.fill_between(x, fit_min, fit_max, color="black", alpha=0.05)

# --- Plot splitting intensity scatter ---
fig3_utils.plot_results_scatter(ax, null_results, "#dddddd", "null", "o")
fig3_utils.plot_results_scatter(ax, average_results, "#777777", "fair", "d")
fig3_utils.plot_results_scatter(ax, good_results, "#333333", "good", "s")

# --- Configure plot details ---
ax.set_ylabel("Splitting intensity")
ax.set_xlim([0, 360])
ax.set_ylim([-1, 1])
ax.yaxis.set_minor_locator(MultipleLocator(0.125))
xtickrange = np.arange(0, 361, 60)
ax.set_xticks(xtickrange)
ax.set_xticklabels([])
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

ax = ax_dict["B"]
fig3_utils.plot_phi_baz(ax, good, "#333333", "Good", "s")
fig3_utils.plot_phi_baz(ax, fair, "#777777", "Fair", "d")

ax.set_ylim([0, 180])
ytickrange = np.arange(0, 181, 45)
ax.set_yticks(ytickrange)
ax.set_yticklabels(ytickrange)
ax.yaxis.set_minor_locator(MultipleLocator(5))
ax.set_ylabel("\u03C6, °")

ax.set_xlim([0, 360])
ax.set_xticks(xtickrange)
ax.set_xticklabels(xtickrange)
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.set_xlabel("Backazimuth, °")
ax.axhline(y=68, linewidth=1, linestyle="--", color="gray")
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

# --- Save figure ---
if not (pathlib.Path.cwd() / "plots").exists():
    (pathlib.Path.cwd() / "plots").mkdir(exist_ok=True, parents=True)
plt.savefig("plots/figure3.pdf", dpi=400, bbox_inches="tight")
