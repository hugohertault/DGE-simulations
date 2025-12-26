#!/usr/bin/env python3
"""
Dark Geometry - Configuration N-body RAMSES/ECOSMOG
====================================================

Ce script génère les fichiers de configuration pour lancer
des simulations N-corps avec Dark Geometry.

Codes supportés :
- RAMSES (avec patch MG)
- ECOSMOG (extension RAMSES pour gravité modifiée)
- MG-Gadget (alternative)
"""

import numpy as np
import os

# =============================================================================
# SIMULATION SPECIFICATIONS
# =============================================================================

SIMULATION_SPECS = {
    # Low resolution (test)
    'low_res': {
        'box_size': 256,       # Mpc/h
        'n_particles': 256,    # 256³
        'n_cells': 256,        # AMR base level
        'levelmax': 14,        # Max refinement
        'z_start': 49,
        'z_end': 0,
        'n_outputs': 10,
    },
    
    # Medium resolution (analysis)
    'med_res': {
        'box_size': 500,
        'n_particles': 512,
        'n_cells': 512,
        'levelmax': 16,
        'z_start': 49,
        'z_end': 0,
        'n_outputs': 20,
    },
    
    # High resolution (production)
    'high_res': {
        'box_size': 1000,
        'n_particles': 1024,
        'n_cells': 1024,
        'levelmax': 18,
        'z_start': 99,
        'z_end': 0,
        'n_outputs': 50,
    },
}

# Cosmological parameters (Planck 2018)
COSMO_PARAMS = {
    'omega_m': 0.315,
    'omega_b': 0.0493,
    'omega_l': 0.685,
    'h': 0.674,
    'sigma8': 0.773,      # DG value!
    'n_s': 0.9649,
}

# Dark Geometry parameters
DG_PARAMS = {
    'alpha_star': 0.075,
    'beta': 2.0/3.0,
    'k_J_eq': 0.05,       # h/Mpc at equality
    'a_eq': 2.93e-4,
    'S_max': 0.882,
}


# =============================================================================
# RAMSES NAMELIST GENERATOR
# =============================================================================

def generate_ramses_namelist(spec_name='med_res', output_dir='./'):
    """
    Generate RAMSES namelist file.
    """
    
    spec = SIMULATION_SPECS[spec_name]
    
    namelist = f"""
&RUN_PARAMS
cosmo=.true.
pic=.true.
poisson=.true.
hydro=.false.
ncontrol=1
nsubcycle=1,1,2,2
verbose=.false.
/

&OUTPUT_PARAMS
foutput={spec['n_outputs']}
aout=0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0
/

&INIT_PARAMS
filetype='grafic'
initfile(1)='./ics'
/

&AMR_PARAMS
levelmin={int(np.log2(spec['n_cells']))}
levelmax={spec['levelmax']}
ngridmax=2000000
npartmax=5000000
boxlen={spec['box_size']:.1f}
/

&POISSON_PARAMS
gravity_type=1
/

&REFINE_PARAMS
m_refine=10*8.
mass_sph=1.d-4
/

&COSMO_PARAMS
omega_m={COSMO_PARAMS['omega_m']}
omega_l={COSMO_PARAMS['omega_l']}
omega_b={COSMO_PARAMS['omega_b']}
h0={COSMO_PARAMS['h']*100}
aexp_ini=0.02
/

! =====================================================
! DARK GEOMETRY PARAMETERS
! =====================================================
&DARK_GEOMETRY_PARAMS
use_dark_geometry=.true.
dg_alpha_star={DG_PARAMS['alpha_star']}
dg_beta={DG_PARAMS['beta']}
dg_k_J_eq={DG_PARAMS['k_J_eq']}
dg_a_eq={DG_PARAMS['a_eq']}
/
"""
    
    filename = os.path.join(output_dir, f'namelist_{spec_name}.nml')
    with open(filename, 'w') as f:
        f.write(namelist)
    
    print(f"Generated: {filename}")
    return filename


# =============================================================================
# ECOSMOG PARAMETER FILE
# =============================================================================

def generate_ecosmog_params(spec_name='med_res', output_dir='./'):
    """
    Generate ECOSMOG parameter file for modified gravity.
    """
    
    spec = SIMULATION_SPECS[spec_name]
    
    params = f"""
# ECOSMOG Parameter File for Dark Geometry
# =========================================

# Simulation type
mg_type = DG           # Dark Geometry

# DG parameters (all derived from first principles)
dg_alpha_star = {DG_PARAMS['alpha_star']}
dg_beta = {DG_PARAMS['beta']}
dg_k_J_eq = {DG_PARAMS['k_J_eq']}
dg_a_eq = {DG_PARAMS['a_eq']}
dg_S_max = {DG_PARAMS['S_max']}

# Solver parameters
poisson_solver = FFT    # or MULTIGRID
fft_type = FFTW
mg_tolerance = 1e-6
mg_max_iter = 100

# Screening (DG doesn't have traditional screening)
screening = none

# Output
output_G_eff = .true.
output_fifth_force = .true.
"""
    
    filename = os.path.join(output_dir, f'ecosmog_dg_{spec_name}.params')
    with open(filename, 'w') as f:
        f.write(params)
    
    print(f"Generated: {filename}")
    return filename


# =============================================================================
# INITIAL CONDITIONS (MUSIC CONFIG)
# =============================================================================

def generate_music_config(spec_name='med_res', output_dir='./'):
    """
    Generate MUSIC configuration for initial conditions.
    """
    
    spec = SIMULATION_SPECS[spec_name]
    
    config = f"""
[setup]
boxlength        = {spec['box_size']}
zstart           = {spec['z_start']}
levelmin         = {int(np.log2(spec['n_particles']))}
levelmax         = {int(np.log2(spec['n_particles']))}
baryons          = no
use_2LPT         = yes
use_LLA          = yes

[cosmology]
Omega_m          = {COSMO_PARAMS['omega_m']}
Omega_L          = {COSMO_PARAMS['omega_l']}
Omega_b          = {COSMO_PARAMS['omega_b']}
H0               = {COSMO_PARAMS['h']*100}
sigma_8          = {COSMO_PARAMS['sigma8']}
nspec            = {COSMO_PARAMS['n_s']}
transfer         = CLASS_DG

[random]
seed[7]          = 12345
seed[8]          = 23456
seed[9]          = 34567

[output]
format           = grafic2
filename         = ./ics

[poisson]
fft_fine         = yes
accuracy         = 1e-5
grad_order       = 6
laplace_order    = 6

# Dark Geometry: Use modified transfer function
[dark_geometry]
use_dg           = yes
alpha_star       = {DG_PARAMS['alpha_star']}
S_max            = {DG_PARAMS['S_max']}
k_J_today        = 0.005
"""
    
    filename = os.path.join(output_dir, f'music_dg_{spec_name}.conf')
    with open(filename, 'w') as f:
        f.write(config)
    
    print(f"Generated: {filename}")
    return filename


# =============================================================================
# SLURM JOB SCRIPT
# =============================================================================

def generate_slurm_script(spec_name='med_res', output_dir='./', 
                          cluster='generic', n_nodes=4):
    """
    Generate SLURM submission script.
    """
    
    spec = SIMULATION_SPECS[spec_name]
    
    # Estimate resources
    if spec_name == 'low_res':
        walltime = '02:00:00'
        n_tasks = 64
    elif spec_name == 'med_res':
        walltime = '24:00:00'
        n_tasks = 256
    else:
        walltime = '72:00:00'
        n_tasks = 1024
    
    script = f"""#!/bin/bash
#SBATCH --job-name=DG_{spec_name}
#SBATCH --output=dg_{spec_name}_%j.out
#SBATCH --error=dg_{spec_name}_%j.err
#SBATCH --nodes={n_nodes}
#SBATCH --ntasks-per-node={n_tasks // n_nodes}
#SBATCH --cpus-per-task=1
#SBATCH --time={walltime}
#SBATCH --partition=standard

# =====================================================
# Dark Geometry N-body Simulation
# Specification: {spec_name}
# Box: {spec['box_size']} Mpc/h
# Particles: {spec['n_particles']}³ = {spec['n_particles']**3:,}
# =====================================================

echo "Starting DG simulation: {spec_name}"
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
mkdir -p output_{spec_name}
cd output_{spec_name}

# Step 1: Generate initial conditions
echo "Generating ICs with MUSIC..."
$HOME/music/MUSIC ../music_dg_{spec_name}.conf

# Step 2: Run RAMSES with DG
echo "Running RAMSES-DG..."
mpirun -np {n_tasks} $RAMSES_DIR/bin/ramses3d ../namelist_{spec_name}.nml

# Step 3: Post-processing
echo "Running halo finder..."
$HOME/rockstar/rockstar -c ../rockstar.cfg

# Step 4: Analysis
echo "Computing power spectrum..."
python ../compute_pk.py --snapshot output_00010/

echo "Done!"
echo "End time: $(date)"
"""
    
    filename = os.path.join(output_dir, f'submit_dg_{spec_name}.sh')
    with open(filename, 'w') as f:
        f.write(script)
    
    os.chmod(filename, 0o755)
    print(f"Generated: {filename}")
    return filename


# =============================================================================
# POWER SPECTRUM COMPUTATION SCRIPT
# =============================================================================

def generate_pk_script(output_dir='./'):
    """
    Generate power spectrum computation script.
    """
    
    script = '''#!/usr/bin/env python3
"""
Compute matter power spectrum from N-body snapshot.
"""

import numpy as np
import argparse
from scipy.fft import fftn
import h5py

def read_snapshot(filename):
    """Read particle positions from snapshot."""
    # Adapt to your format (HDF5, Gadget, RAMSES...)
    with h5py.File(filename, 'r') as f:
        pos = f['PartType1/Coordinates'][:]
        boxsize = f.attrs['BoxSize']
    return pos, boxsize

def compute_density_field(pos, boxsize, ngrid=512):
    """Compute density field on grid using CIC."""
    delta = np.zeros((ngrid, ngrid, ngrid))
    
    # CIC assignment
    cell_size = boxsize / ngrid
    
    for p in pos:
        i = int(p[0] / cell_size) % ngrid
        j = int(p[1] / cell_size) % ngrid
        k = int(p[2] / cell_size) % ngrid
        
        delta[i, j, k] += 1
    
    # Normalize
    mean_density = len(pos) / ngrid**3
    delta = delta / mean_density - 1
    
    return delta

def compute_power_spectrum(delta, boxsize, ngrid):
    """Compute P(k) from density field."""
    # FFT
    delta_k = fftn(delta)
    
    # Power spectrum |δ(k)|²
    pk_3d = np.abs(delta_k)**2 * (boxsize / ngrid)**3
    
    # Bin in k
    kmax = ngrid * np.pi / boxsize
    k_bins = np.linspace(0, kmax, 50)
    k_centers = 0.5 * (k_bins[1:] + k_bins[:-1])
    
    # Compute k for each cell
    kx = np.fft.fftfreq(ngrid, d=boxsize/ngrid) * 2 * np.pi
    ky = np.fft.fftfreq(ngrid, d=boxsize/ngrid) * 2 * np.pi
    kz = np.fft.fftfreq(ngrid, d=boxsize/ngrid) * 2 * np.pi
    
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    K = np.sqrt(KX**2 + KY**2 + KZ**2)
    
    # Bin
    pk = np.zeros(len(k_centers))
    counts = np.zeros(len(k_centers))
    
    for i, (k_lo, k_hi) in enumerate(zip(k_bins[:-1], k_bins[1:])):
        mask = (K >= k_lo) & (K < k_hi)
        if np.sum(mask) > 0:
            pk[i] = np.mean(pk_3d[mask])
            counts[i] = np.sum(mask)
    
    return k_centers, pk

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--snapshot', required=True)
    parser.add_argument('--ngrid', type=int, default=512)
    parser.add_argument('--output', default='pk.txt')
    args = parser.parse_args()
    
    print(f"Reading snapshot: {args.snapshot}")
    pos, boxsize = read_snapshot(args.snapshot)
    
    print(f"Computing density field (ngrid={args.ngrid})...")
    delta = compute_density_field(pos, boxsize, args.ngrid)
    
    print("Computing power spectrum...")
    k, pk = compute_power_spectrum(delta, boxsize, args.ngrid)
    
    # Save
    np.savetxt(args.output, np.column_stack([k, pk]), 
               header='k [h/Mpc]  P(k) [(Mpc/h)³]')
    print(f"Saved: {args.output}")

if __name__ == "__main__":
    main()
'''
    
    filename = os.path.join(output_dir, 'compute_pk.py')
    with open(filename, 'w') as f:
        f.write(script)
    
    os.chmod(filename, 0o755)
    print(f"Generated: {filename}")
    return filename


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    
    print("="*70)
    print("DARK GEOMETRY - N-BODY CONFIGURATION GENERATOR")
    print("="*70)
    
    output_dir = './'
    
    # Generate for all resolutions
    for spec_name in ['low_res', 'med_res', 'high_res']:
        print(f"\n--- {spec_name.upper()} ---")
        spec = SIMULATION_SPECS[spec_name]
        print(f"  Box: {spec['box_size']} Mpc/h")
        print(f"  Particles: {spec['n_particles']}³ = {spec['n_particles']**3:,}")
        print(f"  Memory: ~{spec['n_particles']**3 * 32 / 1e9:.1f} GB (particles only)")
        
        generate_ramses_namelist(spec_name, output_dir)
        generate_ecosmog_params(spec_name, output_dir)
        generate_music_config(spec_name, output_dir)
        generate_slurm_script(spec_name, output_dir)
    
    # Generate common scripts
    generate_pk_script(output_dir)
    
    print("\n" + "="*70)
    print("CONFIGURATION FILES GENERATED")
    print("="*70)
    print("""
To run simulations:

1. Compile RAMSES with DG patch:
   cd ramses/
   patch -p1 < dg_patch.diff
   make

2. Generate initial conditions:
   ./MUSIC music_dg_med_res.conf

3. Run simulation:
   sbatch submit_dg_med_res.sh
   
4. Analyze results:
   python compute_pk.py --snapshot output_00010/
""")
