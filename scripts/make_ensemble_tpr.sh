#!/bin/bash
#SBATCH --job-name=BRER_tpr_setup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=mlerch
#SBATCH --time=01:00:00
#SBATCH --mail-user=jgrosskopf@mcw.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

module load gromacs/2019.6

gmx grompp -f step7_production.mdp -o brer_frame1236_4yay.tpr -c frame_1236_AT1R_4yay_250ns.gro -p topol.top -n index.ndx
