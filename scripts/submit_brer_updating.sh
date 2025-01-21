#!/bin/bash

# example bash script for submitting to a computer cluster using slurm

#SBATCH --job-name=brer_gpu_AT1R_comp4
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4GB
#SBATCH --gres=gpu:1
#SBATCH --account=mlerch
#SBATCH --time=07-00:00:00
#SBATCH --mail-user=jgrosskopf@mcw.edu
#SBATCH --mail-type=ALL
#SBATCH --array=9-16%8

cd $SLURM_SUBMIT_DIR # specify working directory

# load required modules
module load brer
module load gromacs/2019.6
module load miniconda3  #anaconda module

conda activate chilife_env    #anaconda enironment with chiLife

dir=$SLURM_SUBMIT_DIR/ensemble_dir/mem_${SLURM_ARRAY_TASK_ID} #create directory in ensemble_dir directory for each ensemble
tpr_file=$SLURM_SUBMIT_DIR/tpr_files/tpr_files_${SLURM_ARRAY_TASK_ID} #create TPR file collection for each ensemble
mkdir $dir
mkdir $tpr_file

python initialize_updating_protocol.py \
                    --starting-structure frame_1236_AT1R_4yay_250ns.gro \
                    --run-index ${SLURM_ARRAY_TASK_ID}
                    
run_brer run.py \
        --input=$SLURM_SUBMIT_DIR/brer_frame1236_4yay.tpr \
        --workdir=$SLURM_SUBMIT_DIR/ensemble_dir \
        --pairs=$SLURM_SUBMIT_DIR/pair_data/pair_data_${SLURM_ARRAY_TASK_ID}.json \
        --threads-per-sim=$SLURM_CPUS_PER_TASK \
        --ensemble-number=${SLURM_ARRAY_TASK_ID}
        
python run_updating.py \
        --latest-structure $dir \
        --run-index ${SLURM_ARRAY_TASK_ID}

for i in {1..30} #this will specify the number of iterations to do in each ensemble
do
    let rep=i-1
    file=$(find $dir/$rep/production -name "*.gro")
    echo $file
    
    gmx grompp -f step7_production.mdp -o $tpr_file/latest_structure_$rep.tpr -c $file -p topol.top -n index.ndx
    
    run_brer run.py \
            --input=$tpr_file/latest_structure_$rep.tpr \
            --workdir=$SLURM_SUBMIT_DIR/ensemble_dir \
            --pairs=$SLURM_SUBMIT_DIR/pair_data/pair_data_${SLURM_ARRAY_TASK_ID}.json \
            --threads-per-sim=$SLURM_CPUS_PER_TASK \
            --ensemble-number=${SLURM_ARRAY_TASK_ID}
    
    python run_updating.py \
                    --latest-structure $dir \
                    --run-index ${SLURM_ARRAY_TASK_ID}
     
done

