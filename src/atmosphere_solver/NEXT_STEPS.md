# Atmosphere Solver Next Steps

## 1. Finish self-consistent initialization

Complete the current run setup so the initial atmosphere is internally consistent:

- build `z`, `T(z)`, guessed `rho(z)`, `g(z)`, and hydrostatic `P(z)`;
- run EasyChem over the full profile to get species mole fractions;
- compute mean molecular mass and species number densities;
- iterate or refresh `rho`, `P`, and composition until the ideal-gas relation and hydrostatic equation are consistent enough for an initial state.

## 2. Add radiative-convective temperature dynamics

Replace the hardcoded isothermal temperature profile with a first radiative-convective estimate:

- estimate Bond albedo from wavelength-dependent absorption and scattering optical depths;
- compute heating/cooling from stellar flux and atmospheric opacity;
- include convective adjustment when the lapse rate is unstable;
- feed the resulting `T(z)` back into chemistry and hydrostatic structure.

## 3. Split layered atmosphere solvers

Separate the current atmospheric dynamics class into two modes:

- `layered_equilibrium`: local equilibrium chemistry at each layer, useful for initialization and fast tests;
- `layered_disequilibrium`: time-dependent species evolution with transport, photochemistry, escape, and source/sink terms.

## 4. Add planetary-class abundance presets

Add hardcoded initial atomic abundance presets for broad planetary classes:

- Earth-like terrestrial planets;
- Super-Earths;
- mini-Neptunes.

Titan does not need this path for now because the relevant molecular abundances are already known/prescribed.
