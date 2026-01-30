# ProGuide module for modeling proteins from DEER distance distributions

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
    - File name should have the residue pair specified somewhere in it like `_res1_res2_` (ex. `b2ar_148_266_iso_nb80.txt`)
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

3. Within the `initialize_updating_protocol.py` and `run_updating.py` scripts, specify: 
    - The label pairs:
        - `label_pairs = ['148_266']` will specify a single label pair constraint between residues 148 and 266
        - `label_pairs = ['148_266', '148_228']` will specify label pair constraints between 148 and 266, and another set of constraints between 148 and 228. This can be done for any number of label pair constraints. 
    - An identifier for the distance distribution file:
        - If this is a modeling run with a single label pair constraint and the distance distribution file is named `b2ar_148_266_final.txt`: `exp_data = 'b2ar_148_266_final.txt'`
        - if this is a modeling run with multiple label pair constraints and the distribution files are named `b2ar_148_266_final.txt` and `b2ar_148_228_final.txt`, then get an identifier of both and specify like this: `exp_data = glob.glob('*final.txt')`. This will gather all files with `final.txt` in the name. In this case, only the distance distribution files have this and the program will sort which residue pairs to apply the constraints to based on the name of the file.
    - The learning rate:
        - `learn_rate = 0.2` -- this dictates the size of structural change to implement on each model iteration, as a fraction of the calculated correction necessary to overlap the means of the simulated and experimental distribution. I.e. 10$\AA$ is calculated, 10*0.2 = 2$\AA$ structural change implemented in the modeling iteration. 
    - The momentum parameter:
        - `momentum = 1` -- this dictates how many previous model's simulated distribution to take into account when calcluating the structural change/update. A momentum parameter of 1 will calculate the update based only on the model latest model generated. A momentum parameter of 2 will calculate the update based on the last 2 models generated, and so on. This is useful when the experimental distribution is relatively broad compared to the simulated distribution. 
    - **This redundancy to specify these in both files will be fixed in future versions**. Thank you for your understanding :)


4. Specify `A`, `tau`, `tolerance` and `production_time` in `run.py`.
    - `A` -- specifies a searchable range for the bias (see refs: ###). Recommendation: 50-200
    - `tau` -- specifies the amount of time (ps) the system is exposed to the `A` parameter in training before an updated `A` is calculated. Highly depends on system. `tau` = 1000 works well in application tested here.
    - `tolerance` -- specifies the tolerance from the target distance in the convergence phase (in nm). Once the distance is within the tolerance, convergence ends and moves to production.
    -`production_time` -- specifies the amount of time the system relaxes around the bias potential (in ns).

5. Add the starting modeling GRO file, starting model TPR file name, number of replicates, and number of iterations to the `submit_brer_updating.sh` file:
    - This script is written using the SLURM workload manager and using bash. Please click on the file to see a full example of the script to submit and run this modeling.
    - To specify, for example, 8 iterations: `#SBATCH --array=1-8%8`
        - This creates an array job of 8 runs with the specified hardware and modeling parameters. Each replicate will be assigned a folder, and each iteration within each replicate will be assigned a folder. The replicates are numbered by their array number. 
    - Specify the starting structure used to make the TPR file:
        ```
        python initialize_updating_protocol.py \
                    --starting-structure frame_1236_AT1R_4yay_250ns.gro \
                    --run-index ${SLURM_ARRAY_TASK_ID}
        ```
    - Specify the TPR file:
        ```
        run_brer run.py \
            --input=$SLURM_SUBMIT_DIR/brer_frame1236_4yay.tpr \
        ```
    - **You are now ready to run the modeling!**

6. Useful outputs

this is created by constructing a GROMACS MD simulation (can use CHARMM-GUI), (preferably) energy minimizing/equilibrating the system, selecting a frame from the simulation as a .GRO file, and converting to TPR file. The needed files and code to do this is in scripts/make_ensemble_tpr.sh

