# BRER-chiLife module for modeling proteins from DEER distance distributions

## Description
This is open-source code for modeling/refining protein conformations from distance distributions. The examples shown here specifically use distance distributions derived from double electron-electron resonance (DEER) spectroscopy.

## Software Requirements
- [GROMACS 2019.6](https://manual.gromacs.org/documentation/)
- [BRER Plugin](https://github.com/kassonlab/brer_plugin)
- [Run BRER Module](https://github.com/kassonlab/run_brer)
- Python 3.10+
- [chiLife](https://github.com/StollLab/chiLife)

## Usage
### Requirements to get started
- A starting structure (crystal, NMR, CryoEM structure) or model (AlphaFold, Rosetta, etc.)
- Distance distributions between any two spin labels on a protein/protein complex (Note: this can be adapted for other labels or for direct measurements to other atoms in the protein)
- GROMACS TPR file of starting structure -- This involves setting up the starting structure/model in a Gromacs simulation environment ([CHARMM-GUI](https://charmm-gui.org/) can help with this)



this is created by constructing a GROMACS MD simulation (can use CHARMM-GUI), (preferably) energy minimizing/equilibrating the system, selecting a frame from the simulation as a .GRO file, and converting to TPR file. The needed files and code to do this is in scripts/make_ensemble_tpr.sh

