# IAGA Summer School 2019
Materials for the workshop on magnetic observatories and modelling. Tutorials are written in Python as Jupyter Notebooks and can be found in the ``notebooks`` directory.

Authors: David Kerridge, Grace Cox, Ashley Smith (for info, contact ashley.smith@ed.ac.uk)

# Contents
 1. Visualising Geomagnetic Observatory Data
     - plotting observatory time series data, from years to seconds
 2. K Index Calculation
     - how to calculate the activity "K" index
 3. Magnetic Field Line Tracing
     - moving along field lines and finding magnetically conjugate points
 4. Spherical Harmonic Models 1
     - representation and evaluation of the IGRF
 5. Spherical Harmonic Models 2
     - building SH models, and using a satellite-derived virtual observatory (VO) dataset
 6. Accessing Swarm Data and Models Using VirES (viresclient)

 
[Use nbviewer to view the notebooks non-interactively](https://nbviewer.jupyter.org/github/smithara/IAGA_SummerSchool2019/tree/master/notebooks/)

# Running on the ESA Virtual Research Environment (VRE)
The VRE provides cloud-based [JupyterLab](https://jupyterlab.readthedocs.io/) servers where the tutorials can be run interactively. Students at the summer school will be given access.
### To get started:
- Login at https://vre.vires.services/
- Open a terminal (File/New/Terminal)
- (You should be in your home directory)
- `git clone https://github.com/smithara/IAGA_SummerSchool2019.git`

### If a newer version is available:
This will get the newer version and overwrite any local changes:
- `cd ~/IAGA_SummerSchool2019/`
- `git fetch`
- `git reset --hard origin/master`

# Instructions to set up on your own machine
Prerequisite: You will need an installation of Python 3.6+ with recent versions of numpy, matplotlib, pandas, scipy, jupyter. [Anaconda](https://www.anaconda.com/distribution/) is the recommended way to install all this.

Download the project to the location of your choice (the green button near the top right of this page). Or using git, just:
```
git clone https://github.com/smithara/IAGA_SummerSchool2019.git
```

# Project structure

```
.
├── LICENSE
├── README.md
├── data
│   ├── external       <- Data that is sourced from elsewhere
│   ├── downloaded     <- (currently empty)
├── notebooks          <- Jupyter notebooks (tutorials)
└── src                <- Source code for this project (modules)
```

# Included data
- ``data/external/pulsations/``
    - 1Hz observatory data from Hartland (HAD)
    - (two days on which pulsation activity was observed: 14-12-2018 and 11-06-2019)
    - figures showing the pulsations on these days. 
- ``data/external/k_inds/``
    - Observatory K indices
    - http://geomag.bgs.ac.uk/data_service/data/magnetic_indices/k_indices.html
- ``data/external/ESK_hourly/``
- ``data/external/ESK_2003/``
    - INTERMAGNET 1-minute data from Eskdalemuir observatory in IAGA-2002 format
    - http://www.intermagnet.org/data-donnee/download-eng.php
- ``data/external/SwarmVO_IAGASummerSchool.dat
    - Virtual Observatory (VO) dataset built from Swarm data
- ``data/external/igrf12coeffs.txt``
    - [International Geomagnetic Reference Field 12 (IGRF-12)](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
- ``data/external/land_5deg.csv``
    - A simple binary map of the world with 5-degree resolution
- ``data/external/oamjumpsapplied.dat``
-  ``data/external/quiet_days.txt``


# Useful links
- MagPySV:
    - https://magpysv.readthedocs.io
    - https://github.com/gracecox/MagPySV-examples
    - https://doi.org/10.1029/2018GC007714
- viresclient:
    - http://viresclient.readthedocs.io

- Setting up a structured repository for reproducible science:
    - https://github.com/mkrapp/cookiecutter-reproducible-science