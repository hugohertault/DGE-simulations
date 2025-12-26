#!/bin/bash
#SBATCH --job-name=DG_low_res
#SBATCH --output=dg_low_res_%j.out
#SBATCH --error=dg_low_res_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --partition=standard

# =====================================================
# Dark Geometry N-body Simulation
# Specification: low_res
# Box: 256 Mpc/h
# Particles: 256³ = 16,777,216
# =====================================================

echo "Starting DG simulation: low_res"
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
mkdir -p output_low_res
cd output_low_res

# Step 1: Generate initial conditions
echo "Generating ICs with MUSIC..."
$HOME/music/MUSIC ../music_dg_low_res.conf

# Step 2: Run RAMSES with DG
echo "Running RAMSES-DG..."
mpirun -np 64 $RAMSES_DIR/bin/ramses3d ../namelist_low_res.nml

# Step 3: Post-processing
echo "Running halo finder..."
$HOME/rockstar/rockstar -c ../rockstar.cfg

# Step 4: Analysis
echo "Computing power spectrum..."
python ../compute_pk.py --snapshot output_00010/

echo "Done!"
echo "End time: $(date)"
