# Supplementary analysis/visualisation code
Data and code required to reproduce the visualisations presented in the manuscript:

Lithospheric anisotropy characterises the post-subduction setting of northern Borneo: new results from XKS splitting analysis

Bacon, C.A., Rawlinson, N., Pilia, S., Gilligan, A., Wehner, D., Cornwell, D.G., Tongkul, F.

Pre-print: [![DOI](https://img.shields.io/badge/EarthArxiv-10.31223/X57H14-blue)](https://doi.org/10.31223/X57H14)

Supplementary datafiles: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6461787.svg)](https://doi.org/10.5281/zenodo.6461787)

Analysis and visualisation: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6480581.svg)](https://doi.org/10.5281/zenodo.6480581)

## Steps to reproduce figures
1. Clone this repository and navigate to it, e.g.:

```
git clone https://github.com/hemmelig/g-cubed-manuscript
cd g-cubed-manuscript
```

2. Install the packages listed in the environment.yml file, either manually, or using (for example) conda:

```
conda env create
```

3. Add `splitracer_summary.mplstyle` (a matplotlib stylesheet) to your `mpl_configdir` (usually found at `~/.config/matplotlib`)

4. Optional: Install Helvetica font for Matplotlib

5. Run `download_grd_datafiles.gmt` script to download the DEM data used in Figure 1:

```
cd figures
bash download_grd_datafiles.gmt
```

6. Navigate to each figure directory and run the `.gmt` (as `bash <script>.gmt`) or `.py` (as `python <script>.py`) scripts.

## Notes
These figures were prepared using Linux 20.04. A limited number were produced using Affinity Designer for iPad, a licensed piece of software for graphic design. `.afdesign` files are provided.

For Supplementary Figure S6, you will need `AnisotroPy`, an open-source package for the analysis of seismic anisotropy. This can be downloaded from GitHub at https://github.com/hemmelig/AnisotroPy, where you will also find instructions on how to install the package.
