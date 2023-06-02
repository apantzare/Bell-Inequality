#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=A100:4
#SBATCH -A SNIC2022-5-398
#SBATCH -t 15:00:00
#SBATCH --job-name=Axels_test

module purge
LMOD_DISABLE_SAME_NAME_AUTOSWAP='no'
module load scikit-learn
module load TensorFlow/2.11
#module load SciPy-bundle/2020.11-fosscuda-2020b

# Copy data to TMPDIR
cd $TMPDIR
cp /cephyr/users/paaxel/Alvis/SVC-XGBoost.py .
cp /cephyr/users/paaxel/Alvis/data.dat .
cp /cephyr/users/paaxel/Alvis/target.dat .
cp /cephyr/users/paaxel/Alvis/Data.csv .

# Run your code
python3 SVC-XGBoost.py
# Retrieve results
cp result_SVC.dat $SLURM_SUBMIT_DIR
cp pred_SVC.dat $SLURM_SUBMIT_DIR
cp result_real_SVC.dat $SLURM_SUBMIT_DIR

cp SVC.pickle $SLURM_SUBMIT_DIR
