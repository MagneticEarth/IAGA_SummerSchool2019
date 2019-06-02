# IAGA Summer School 2019
Materials for the workshop on magnetic observatories and modelling. Tutorials are written in Python as Jupyter Notebooks and can be found in the ``notebooks`` directory:
 1. (tutorial name)
    - (description)
 2. 
 ...
 
[Use nbviewer to view the notebooks non-interactively](https://nbviewer.jupyter.org/github/smithara/IAGA_SummerSchool2019/tree/master/notebooks/)

# Running on the ESA Virtual Research Environment (VRE)
The VRE provides cloud-based [JupyterLab](https://jupyterlab.readthedocs.io/en/stable/]) servers where the tutorials can be run interactively. Students at the summer school will be given access.
- [VirES for Swarm](https://vires.services/)
- [VRE](https://vre.vires.services/)
...

# Instructions to set up on your own machine
Prerequisite: You will need an installation of Python 3.6+ with recent versions of numpy, matplotlib, pandas, scipy, jupyter. (link to instructions to set this up using conda)

Download the project to the location of your choice. With git, just:
```
git clone https://github.com/smithara/IAGA_SummerSchool2019.git
```
...


# Data sources
 - ``data/external/igrf12coeffs.txt``
     - [International Geomagnetic Reference Field 12 (IGRF-12)](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
 - ``data/external/ESK_2003/``  [not in repository]
     - INTERMAGNET 1-minute data from Eskdalemuir observatory in IAGA-2002 format
     - I tried to download this myself from http://www.intermagnet.org/data-donnee/download-eng.php but I just got an error
 - ``data/external/k_inds/``  [not in repository]
     - Observatory K indices
     - http://geomag.bgs.ac.uk/data_service/data/magnetic_indices/k_indices.html ?
 - ``data/external/land_5deg.csv``
     - A simple binary map of the world with 5-degree resolution

# Project structure

```
.
├── LICENSE
├── README.md
├── data
│   ├── external       <- Data that is sourced from elsewhere (these can be part of the repository if light)
│   ├── downloaded     <- Data that can be automatically downloaded by a script within src
├── notebooks          <- Jupyter notebooks (tutorials)
└── src                <- Source code for this project (modules)
```