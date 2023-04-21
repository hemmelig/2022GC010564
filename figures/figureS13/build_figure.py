# -*- coding: utf-8 -*-
"""
This script can be used to generate supplemental Figure S13 of:

    Bacon, C. A., Rawlinson, N., Pilia, S., Gilligan, A., Wehner, D.,
    Cornwell, D. G., & Tongkul, F. (2022). The signature of lithospheric
    anisotropy at post-subduction continental margins: New insight from
    XKS splitting analysis in northern Borneo. Geochemistry, Geophysics,
    Geosystems, 23, e2022GC010564. https://doi.org/10.1029/2022GC010564

"""

import subprocess

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import pandas as pd


plt.style.use("2022GC010564")
mpl.rcParams["font.family"] = "Helvetica"

cmd = (
    "awk '{print $2,$3}' "
    "../data/splitting_results/supp_file_7_1layer_averages.tab | gmt grdtrack "
    "-G../gmt_datafiles/grd/thickness_LAB.grd > lab_thickness_at_stations.txt"
)
subprocess.run(cmd, shell=True)

df1 = pd.read_csv("lab_thickness_at_stations.txt", delimiter="\t", header=None)
df2 = pd.read_csv(
    "../data/splitting_results/supp_file_7_1layer_averages.tab",
    delim_whitespace=True,
    header=None,
)

df2[5] = df1[2]
df2 = df2.dropna()

fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(df2[5].values, df2[4].values, marker="+", c="k")
ax.set_xlabel("LAB depth, km")
ax.set_ylabel("Delay time, s")

corr, _ = pearsonr(df2[5].values, df2[4].values)

ax.text(
    0.2,
    0.9,
    f"Pearson R = {corr:5.3f}",
    fontsize=16,
    va="center",
    ha="center",
    transform=ax.transAxes,
)

plt.savefig("figureS13.pdf", bbox_inches="tight")
