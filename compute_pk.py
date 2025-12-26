#!/usr/bin/env python3
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
