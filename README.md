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
- toppar folder with ITP files (output from CHARMM-GUI)

### A simple single constraint example
This is a simple example for the use of a single constraint (single distance distribution), though a muli-constraint (multi-distribution) example is very similar! The below instructions assume that all files mentioned are within the same folder/directory.
1. Create the TPR file of the starting structure simulation environment (Note: it is recommended that minimization and equilibration is done prior to creating the TPR file).
    - In the make_ensemble_tpr.sh file, you will specify the following:
    ```
    gmx grompp -f simulation.mdp -o output_tpr_file.tpr -c input_for_tpr.gro -p topol.top -n index.ndx
    ```
    where 
    `-f` species the simulation instructions. Make sure your simulation time exceeds production time x iterations.
    `-o` is the name you choose for the TPR file.
    `-c` is the input file. This must be a `.gro` file from the previous minimization/equilibration. 
    `-p` is the name of the topology list file, specifies topology files from toppar folder
    `-n` is the index file for the input file, species atom indices

2. Add distance distribution data file to your working directory.
    - File name should have the residue pair specified somewhere in it (ex. `b2ar_148_266_iso_nb80.txt`)
    - Data in file should be formatted in columns where the first column is distance (in angstroms) and  the second column is the associated probability, with **no labels/data headers**:
        ```
        1.505000000000000071e+01 4.552121762066039048e-28
        1.555000000000000071e+01 2.541999450004171662e-27
        1.605000000000000071e+01 1.384918043307330848e-26
        ...
        7.904999999999999716e+01 1.564411647338884091e-19
        7.954999999999999716e+01 3.715811214973827870e-20
        8.004999999999999716e+01 8.610795612490454540e-21
        ```

this is created by constructing a GROMACS MD simulation (can use CHARMM-GUI), (preferably) energy minimizing/equilibrating the system, selecting a frame from the simulation as a .GRO file, and converting to TPR file. The needed files and code to do this is in scripts/make_ensemble_tpr.sh

