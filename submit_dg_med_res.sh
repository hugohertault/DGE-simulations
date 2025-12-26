#!/bin/bash
#SBATCH --job-name=DG_med_res
#SBATCH --output=dg_med_res_%j.out
#SBATCH --error=dg_med_res_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --partition=standard

# =====================================================
# Dark Geometry N-body Simulation
# Specification: med_res
# Box: 500 Mpc/h
# Particles: 512³ = 134,217,728
# =====================================================

echo "Starting DG simulation: med_res"
echo "Date: $(date)"
echo "Host: $(hostname)"

# Load modules
module load intel/2021
module load mpi/openmpi-4.1
module load fftw/3.3.10

# Set up environment
export OMP_NUM_THREADS=1
export RAMSES_DIR=$HOME/ramses-dg

# Create output directory
mkdir -p output_med_res
cd output_med_res

# Step 1: Generate initial conditions
echo "Generating ICs with MUSIC..."
$HOME/music/MUSIC ../music_dg_med_res.conf

# Step 2: Run RAMSES with DG
echo "Running RAMSES-DG..."
mpirun -np 256 $RAMSES_DIR/bin/ramses3d ../namelist_med_res.nml

# Step 3: Post-processing
echo "Running halo finder..."
$HOME/rockstar/rockstar -c ../rockstar.cfg

# Step 4: Analysis
echo "Computing power spectrum..."
python ../compute_pk.py --snapshot output_00010/

echo "Done!"
echo "End time: $(date)"
