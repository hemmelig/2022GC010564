# -*- coding: utf-8 -*-
"""
This script will create a summary figure for an example splitting measurement
Specifically, this script can be used to generate supplemental Figure S4 of:

    Lithospheric anisotropy characterises the post-subduction setting of
    northern Borneo: new results from XKS splitting analysis
    Bacon, C.A., Rawlinson, N., Pilia, S., Gilligan, A., Wehner, D.,
    Cornwell, D.G., Tongkul, F.

which has been submitted to Earth and Planetary Science Letters.

"""

import pathlib

import functions as splitracer_tools


results_dir = pathlib.Path.cwd() / "data"
event = "49"

splitracer_tools.plot_summary(results_dir, event, phi=71.0, dt=0.80, null=False)
