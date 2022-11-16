# -*- coding: utf-8 -*-
"""
This script will create a summary figure for an example splitting measurement
Specifically, this script can be used to generate supplemental Figure S5 of:

    Bacon, C. A., Rawlinson, N., Pilia, S., Gilligan, A., Wehner, D.,
    Cornwell, D. G., & Tongkul, F. (2022). The signature of lithospheric
    anisotropy at post-subduction continental margins: New insight from
    XKS splitting analysis in northern Borneo. Geochemistry, Geophysics,
    Geosystems, 23, e2022GC010564. https://doi.org/10.1029/2022GC010564

"""

import pathlib

import functions as splitracer_tools


results_dir = pathlib.Path.cwd() / "data"
event = "18"

splitracer_tools.plot_summary(results_dir, event, phi=56, dt=1.03, null=True)
