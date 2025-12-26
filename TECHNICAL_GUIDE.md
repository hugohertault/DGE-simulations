# Dark Geometry CLASS - Guide Technique

## Installation

### Prérequis

```bash
# Compilateur C/C++
sudo apt install gcc g++ make

# Python et dépendances
pip install numpy scipy matplotlib cython emcee corner
```

### Compilation

```bash
# Extraire l'archive
tar -xzf class_dg_complete.tar.gz
cd class_dg

# Compiler CLASS
make clean
make class

# Compiler le binding Python
cd python
python setup.py build_ext --inplace
cd ..

# Créer les liens symboliques pour les données
mkdir -p python/external
ln -sf $(pwd)/external/bbn python/external/bbn
ln -sf $(pwd)/external/HyRec2020 python/external/HyRec2020
```

---

## Utilisation en Ligne de Commande

### Exécuter CLASS-DG

```bash
# Avec Dark Geometry activé
./class dg_test.ini

# Référence ΛCDM (pour comparaison)
./class lcdm_ref.ini
```

### Fichiers de Sortie

```
output/
├── dg_test_00_pk.dat      # P(k) à z=0
├── dg_test_00_pk_cb.dat   # P(k) CDM+baryons
├── dg_test_00_cl.dat      # C_ℓ CMB
└── dg_test_00_cl_lensed.dat
```

---

## Utilisation Python

### Calcul Simple

```python
import sys
sys.path.insert(0, '/path/to/class_dg/python')
import os
os.chdir('/path/to/class_dg')

from classy import Class

# Initialiser
cosmo = Class()
cosmo.set({
    'output': 'mPk,tCl,pCl,lCl',
    'lensing': 'yes',
    'P_k_max_h/Mpc': 10,
    'l_max_scalars': 2500,
    'omega_b': 0.02237,
    'omega_cdm': 0.1200,
    'H0': 67.36,
    'A_s': 2.1e-9,
    'n_s': 0.9649,
    'tau_reio': 0.0544,
})

# Calculer
cosmo.compute()

# Extraire les résultats
sigma8 = cosmo.sigma8()
pk = cosmo.pk(0.1, 0)  # P(k=0.1 h/Mpc, z=0)
cls = cosmo.lensed_cl(2500)

print(f"σ₈ = {sigma8:.4f}")
print(f"P(k=0.1) = {pk:.2e}")

# Nettoyer
cosmo.struct_cleanup()
cosmo.empty()
```

### Calculer le Spectre Complet

```python
import numpy as np

k_values = np.logspace(-4, 1, 200)  # h/Mpc
pk_values = np.array([cosmo.pk(k, 0) for k in k_values])

# Sauvegarder
np.savetxt('pk_dg.dat', np.column_stack([k_values, pk_values]),
           header='k [h/Mpc]  P(k) [(Mpc/h)^3]')
```

### Calculer σ₈ Manuellement

```python
def compute_sigma8(k, pk, R=8.0):
    """Calcul de σ₈ par intégration."""
    def W(x):
        # Filtre top-hat
        return np.where(np.abs(x) > 1e-4,
                       3*(np.sin(x) - x*np.cos(x))/x**3,
                       1.0)
    
    integrand = k**2 * pk * W(k*R)**2 / (2*np.pi**2)
    return np.sqrt(np.trapz(integrand, k))
```

---

## Structure du Code DG

### Fichier Principal : `dark_geometry.c`

```c
// Constantes fondamentales
static const double S_MAX = 0.882;

// Fonction de suppression
double dg_suppression_factor(double k, double a, 
                             struct dark_geometry_params * pdg) {
    double k_J = dg_k_Jeans(a, pdg);
    double x2 = (k / k_J) * (k / k_J);
    return (1.0 + S_MAX * x2) / (1.0 + x2);
}

// Échelle de Jeans
double dg_k_Jeans(double a, struct dark_geometry_params * pdg) {
    // Évolution avec a...
    return k_J;
}
```

### Application dans `fourier.c`

```c
// Après le calcul de P(k) ΛCDM
if (ppt->dg.has_dg == _TRUE_) {
    for (index_k = 0; index_k < pfo->k_size; index_k++) {
        k = pfo->k[index_k];
        S = dg_suppression_factor(k, a, &(ppt->dg));
        // ln(P_DG) = ln(P_ΛCDM) + ln(S)
        pfo->ln_pk_l[...] += log(S);
    }
}
```

---

## Paramètres DG

### Dans `dark_geometry.h`

```c
#define _DG_ALPHA_STAR_UV_ 0.075   // Couplage UV
#define _DG_ALPHA_STAR_IR_ 0.083   // Couplage IR
#define _DG_BETA_ 0.666666666667   // Exposant holographique
#define _DG_XI_ 0.10               // Couplage non-minimal
#define _DG_K_J_EQ_ 0.05           // k_J à l'égalité [h/Mpc]
```

### Modifier les Paramètres

Pour tester différentes valeurs, éditer `dark_geometry.h` et recompiler :

```bash
# Éditer
vim include/dark_geometry.h

# Recompiler
make clean && make class
```

---

## MCMC

### Lancer le MCMC

```bash
python run_mcmc.py
```

### Configuration

Dans `run_mcmc.py` :

```python
# Paramètres MCMC
nwalkers = 16      # Nombre de walkers (≥ 2×ndim)
nsteps = 500       # Nombre de pas
burnin = nsteps//4 # Burn-in

# Likelihood
SIGMA8_OBS = 0.762  # Valeur observée
SIGMA8_ERR = 0.020  # Erreur
```

### Résultats

```
output/
├── mcmc_results.png   # Histogramme σ₈
└── mcmc_corner.png    # Triangle plot
```

---

## Validation CMB

### Lancer la Validation

```bash
python cmb_validation.py
```

### Sortie

```
output/
└── cmb_validation.png  # Spectres TT, EE, résidus
```

---

## Dépannage

### Erreur : fichiers BBN non trouvés

```bash
# Créer les liens symboliques
ln -sf $(pwd)/external/bbn python/external/bbn
ln -sf $(pwd)/external/HyRec2020 python/external/HyRec2020
```

### Erreur : classy non trouvé

```bash
# Recompiler le binding Python
cd python
python setup.py build_ext --inplace
```

### σ₈ incorrect

Vérifier que `has_dg = _TRUE_` dans `dark_geometry_init()`.

---

## Référence API

### Fonctions DG

| Fonction | Description |
|----------|-------------|
| `dg_suppression_factor(k,a,pdg)` | S(k) ∈ [0.882, 1] |
| `dg_k_Jeans(a,pdg)` | k_J(a) en h/Mpc |
| `dg_alpha_star(k,pdg)` | α*(k) avec running |
| `dg_G_eff_ratio(k,a,pdg)` | G_eff/G ≥ 1 |
| `dg_mu(k,a,pdg)` | μ pour Poisson (=1) |

### Structure `dark_geometry_params`

```c
struct dark_geometry_params {
    short has_dg;           // Activation DG
    short has_dg_e;         // Extension DG-E
    double alpha_star_UV;   // Couplage UV
    double beta;            // Exposant holographique
    double k_J_eq;          // k_J à l'égalité
    double Omega_m_0;       // Densité matière
    double h;               // Hubble réduit
};
```

---

*Guide technique - Dark Geometry CLASS v1.0*
