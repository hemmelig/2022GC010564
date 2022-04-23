# Lithospheric anisotropy characterises the post-subduction setting of northern Borneo: new results from XKS splitting analysis
Data and code required to reproduce the results presented in the manuscript:

Lithospheric anisotropy characterises the post-subduction setting of northern Borneo: new results from XKS splitting analysis

Bacon, C.A., Rawlinson, N., Pilia, S., Gilligan, A., Wehner, D., Cornwell, D.G., Tongkul, F.

## Steps to reproduce figures
1. Install the packages listed in the environment.yml file, either manually, or using (for example) conda:

```
conda install -f environment.yml
```

2. Add `splitracer_summary.mplstyle` (a matplotlib stylesheet) to your `mpl_configdir` (usually found at `~/.config/matplotlib`)

3. Optional: Install Helvetica font for Matplotlib

4. Run `download_grd_datafiles.gmt` script to download the DEM data used in Figure 1:

```
cd figures
bash download_grd_datafiles.gmt
```

5. Navigate to each figure directory and run the `.gmt` or `.py` scripts.

## Notes
These figures were prepared using Linux 20.04. A limited number were produced using Affinity Designer for iPad, a licensed piece of software for graphic design. `.afdesign` files are provided.

For Supplementary Figure S6, one must install `AnisotroPy`, an open-source package for the analysis of seismic anisotropy. This can be downloaded from GitHub at https://github.com/hemmelig/AnisotroPy, where one can also find instructions on how to install the package.
