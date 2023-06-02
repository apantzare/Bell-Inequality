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


# Copy data to TMPDIR
cd $TMPDIR
cp /cephyr/users/paaxel/Alvis/1Dconv.py .
cp /cephyr/users/paaxel/Alvis/data.dat .
cp /cephyr/users/paaxel/Alvis/target.dat .
cp /cephyr/users/paaxel/Alvis/Data.csv .

# Run your code
python3 1Dconv.py
# Retrieve results
cp result1Dconv.dat $SLURM_SUBMIT_DIR
#cp true1Dconv.dat $SLURM_SUBMIT_DIR
cp pred1Dconv.dat $SLURM_SUBMIT_DIR
cp result_real_1Dconv.dat $SLURM_SUBMIT_DIR

cp 1Dconv_128+256+128_00005learning_.hdf5 $SLURM_SUBMIT_DIR
