# -*- coding: utf-8 -*-
"""
Helper functions for generating Python summaries of SplitRacer results.

Author: Conor Bacon
Date: 29 March 2022
"""

import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import obspy
import pandas as pd

try:
    plt.style.use("2022GC010564")
except OSError:
    warnings.warn(
        "You have not added the '2022GC010564.mplstyle' stylesheet to your"
        " matplotlib config - continuing with matplotlib defaults."
    )
mpl.rcParams["font.family"] = "Helvetica"


# --- i/o function space ---
def resolve(parent, name):
    return list(parent.glob(name))[0]


def read_energy_grid(path, event):
    """
    Read in the files containing the energy grid.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).

    Returns
    -------
    energy : numpy.ndarray of float, shape(3, n_phi*n_dt)
        Energy for each (phi, dt) pair.

    """

    return pd.read_csv(
        path / f"event{event}_energy_grid.txt", delim_whitespace=True, header=None
    ).values.T


def read_hodograms(path, event):
    """
    Read in the files containing particle motion information.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).

    Returns
    -------
    original_pm : numpy.ndarray of float, shape(2, n_samples)
        Uncorrected particle motion.
    corrected_pm : numpy.ndarray of float, shape(2, n_samples)
        Corrected particle motion.

    """

    original_pm = pd.read_csv(
        path / f"event{event}_RT_orig_pm.txt", delim_whitespace=True, header=None
    )
    corrected_pm = pd.read_csv(
        path / f"event{event}_RT_corr_pm.txt", delim_whitespace=True, header=None
    )

    return original_pm.values.T, corrected_pm.values.T


def read_waveforms(path, event):
    """
    Read in the waveform files.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the origin time).

    Returns
    -------
    st : obspy.Stream object
       Waveform data in the form of a Stream of Traces.

    """

    times, n, e = pd.read_csv(
        path / f"event{event}_tiNE.txt", delim_whitespace=True, header=None
    ).values.T

    st = obspy.Stream()
    for comp, data in zip("NE", [n, e]):
        tr = obspy.Trace()
        tr.data = np.array(data)
        tr.stats.sampling_rate = 1 / (times[1] - times[0])
        tr.stats.component = comp.upper()
        st += tr

    return st


def read_rotated_waveforms(path, event):
    """
    Read in the rotated waveform files.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the origin time).

    Returns
    -------
    st : `obspy.Stream` object
        Stream containing the radial and transverse components of the seismic
        data.

    """

    times, r, t = pd.read_csv(
        path / f"event{event}_tiRT.txt", delim_whitespace=True, header=None
    ).values.T

    st = obspy.Stream()
    for comp, data in zip("RT", [r, t]):
        tr = obspy.Trace()
        tr.data = np.array(data)
        tr.stats.sampling_rate = 1 / (times[1] - times[0])
        tr.stats.component = comp.upper()
        st += tr

    return st


def plot_summary(path, event, dt=0.0, phi=0.0, null=False):
    """
    Utility function bringing together all of the plotting methods for each
    panel.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).
    dt : float, optional
        Measured optimal delay time.
    phi : float, optional
        Measured optimal fast orientation.
    null : bool, optional
        Toggle to cancel plotting of optimal splitting parameters cross.

    """

    # --- Create figure and add axes ---
    fig = _build_grid()

    # --- Filtered waveforms ---
    st_ne = read_waveforms(path, event)
    _filtered_waveforms(fig.axes[0], st_ne, "a")

    # --- R/T waveforms ---
    st_rt = read_rotated_waveforms(path, event)
    _rotated_waveforms(fig.axes[2], st_rt, "b")

    # --- Hodograms ---
    opm, cpm = read_hodograms(path, event)
    _hodogram(fig.axes[1], opm[0], opm[1], "c")
    _hodogram(fig.axes[3], cpm[0], cpm[1], "d")

    # --- Energy grid ---
    dts, phis, vals = read_energy_grid(path, event)
    _energy_grid(fig.axes[4], dts, phis, vals, "e")

    if not null:
        fig.axes[4].scatter(dt, phi, marker="+", label="Best-fitting (\u03C6, \u03B4t)")

    plt.savefig("figureS5.pdf", dpi=400, bbox_inches="tight")


def _build_grid():
    """
    Utility function to construct the grid of panels in the figure.

    Axes are:
        0: Filtered waveforms
        1: Parallel/Perpendicular waveforms before and after correction
      2-5: hodograms
        6: Energy grid

    Returns
    -------
    fig : matplotlib.pyplot.Figure object
        A fully-specified figure containing a grid of panels.

    """

    fig = plt.figure(figsize=(30, 12))
    grid_specs = {"nrows": 6, "ncols": 15, "wspace": 2, "hspace": 2}

    for i in [0, 3]:
        spec = GridSpec(**grid_specs).new_subplotspec((i, 0), colspan=6, rowspan=3)
        fig.add_subplot(spec)
        spec = GridSpec(**grid_specs).new_subplotspec((i, 6), colspan=3, rowspan=3)
        fig.add_subplot(spec)

    spec = GridSpec(**grid_specs).new_subplotspec((0, 9), colspan=6, rowspan=6)
    fig.add_subplot(spec)

    return fig


def _filtered_waveforms(ax, st, letter):
    """
    Plot the filtered waveforms as ZNE components.

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes on which to plot waveforms.
    st : obspy.Stream object
        Stream containing the filtered waveforms.
    letter : str
        Panel identifier.

    """

    # Taper the data
    for tr in st:
        tr.taper(0.15)

    norm = 1.5 * max([max(abs(tr.data)) for tr in st])
    for i, comp in enumerate("EN"):
        tr = st.select(component=comp)[0]
        ax.plot(tr.times(), tr.data / norm + 2 * i, linewidth=1, zorder=2)

    ax.vlines(50, -1.5, 3.5, linewidth=0.8, ls="--", color="k", alpha=0.6)

    ax.set_xlim([st[0].times()[0], st[0].times()[-1]])
    ax.set_ylim([-1.1, 3.1])
    ax.set_yticks(range(0, 3, 2))
    ax.set_yticklabels("EN")
    ax.set_xlabel("Seconds")
    ax.set_ylabel("Component")
    ax.text(
        0.035,
        0.94,
        f"({letter})",
        va="center",
        ha="center",
        fontweight="bold",
        transform=ax.transAxes,
        fontsize=18,
    )


def _rotated_waveforms(ax, st, letter):
    """
    Plot the filtered waveforms as ZNE components.

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes on which to plot waveforms.
    st : obspy.Stream object
        Stream containing the filtered waveforms.
    letter : str
        Panel identifier.

    """

    # Taper the data
    for tr in st:
        tr.taper(0.15)

    norm = 1.5 * max([max(abs(tr.data)) for tr in st])
    for i, comp in enumerate("TR"):
        tr = st.select(component=comp)[0]
        ax.plot(tr.times(), tr.data / norm + 2 * i, linewidth=1, zorder=2)

    ax.vlines(50, -1.5, 3.5, linewidth=0.8, ls="--", color="k", alpha=0.6)

    ax.set_xlim([st[0].times()[0], st[0].times()[-1]])
    ax.set_ylim([-1.1, 3.1])
    ax.set_yticks(range(0, 3, 2))
    ax.set_yticklabels("TR")
    ax.set_xlabel("Seconds")
    ax.set_ylabel("Component")
    ax.text(
        0.035,
        0.94,
        f"({letter})",
        va="center",
        ha="center",
        fontweight="bold",
        transform=ax.transAxes,
        fontsize=18,
    )


def _hodogram(ax, x, y, letter):
    """
    Plot the fast and slow components.

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes on which to plot fast and slow waveforms.
    x : numpy.ndarray of float, shape(n_samples)
        Component 1 for particle motion diagram.
    y : numpy.ndarray of float, shape(n_samples)
        Component 2 for particle motion diagram.
    letter : string
        Panel identifier.

    """

    ax.plot(x, y, linewidth=1, color="k", alpha=0.9)
    ax.set_xlim([-1, 1])
    ax.set_xlabel("Transverse")
    ax.set_xticklabels([])
    ax.set_ylim([-1, 1])
    ax.set_ylabel("Radial")
    ax.set_yticklabels([])
    ax.text(
        0.07,
        0.94,
        f"({letter})",
        va="center",
        ha="center",
        fontweight="bold",
        transform=ax.transAxes,
        fontsize=18,
    )


def _energy_grid(ax, dts, phis, vals, letter):
    """
    Plot the parameter search grid.

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes on which to plot search grid results.
    dts : numpy.ndarray of float, shape(n_dt)
        Range of dt values over which search was performed.
    phis : numpy.ndarray of float, shape(n_phi)
        Range of phi values over which search was performed.
    vals : numpy.ndarray of float, shape(n_dt * n_phi)
        Results of splitting correction.
    letter : str
        Panel identifier.

    """

    # Parse into X, Y, and Z
    Z = vals.reshape(len(set(phis)), len(set(dts))).T
    X, Y = np.mgrid[0.05 : 3.95 : Z.shape[0] * 1j, -90 : 90 : Z.shape[1] * 1j]

    ax.pcolormesh(X, Y, Z, edgecolors="face", cmap="inferno_r")
    ax.contour(X, Y, Z, 9, colors="#cccccc")

    ax.set_ylabel("Fast orientation, °")
    ax.set_xlabel("Delay time, s")
    ax.set_xlim([0, 4])
    ax.set_ylim([-90, 90])
    ax.set_yticks(range(-90, 91, 15))
    ax.set_yticklabels(range(-90, 91, 15))
    ax.text(
        0.035,
        0.97,
        f"({letter})",
        va="center",
        ha="center",
        fontweight="bold",
        transform=ax.transAxes,
        fontsize=18,
    )
