# Solar System Atmosphere Presets

Goal for next work: define a small set of Solar System benchmark atmosphere models that provide:

- planet or moon mass `M_planet`;
- planet or moon radius `R_planet`;
- atmospheric mass fraction `f_atm = M_atm / M_planet`;
- atomic abundances for the equilibrium chemistry input;
- chemical species for the equilibrium chemistry input.

These presets should be used first to test the layered equilibrium atmosphere solver, then later to define broad exoplanet atmosphere classes.

## Notes

- `f_atm` is the normalization needed if surface pressure is meant to be predicted rather than imposed.
- Atmospheric mass is a derived quantity:

```text
M_atm = f_atm * M_planet
```

- For gas giants and ice giants, `f_atm` is not directly comparable because the atmosphere/envelope is a major part of the object. Use a reference-pressure level instead, or treat them as a separate class.
- Atomic abundances below are approximate starting points. They should be derived from the intended molecular composition and then normalized consistently for EasyChem.

## Candidate Presets

| Body | `M_planet` kg | `R_planet` m | `f_atm` | Dominant molecules | Atomic abundance sketch |
|---|---:|---:|---:|---|---|
| Earth | `5.972e24` | `6.371e6` | `8.6e-7` | N2, O2, H2O, CO2 | N/O dominated; include H, C trace |
| Mars | `6.417e23` | `3.390e6` | `3.9e-8` | CO2, N2, Ar | C/O dominated; include N trace |
| Venus | `4.867e24` | `6.052e6` | `9.9e-5` | CO2, N2 | C/O dominated; include N |
| Titan | `1.345e23` | `2.575e6` | `6.7e-5` | N2, CH4, H2 | N/H/C dominated |
| Jupiter reference | `1.898e27` | `6.991e7` | TBD | H2, He, CH4, NH3, H2O | solar-like H/He with C/N/O traces |
| Neptune reference | `1.024e26` | `2.462e7` | TBD | H2, He, CH4 | H/He dominated with enhanced C |

## Implementation Plan

1. Add a preset module, for example `src/planet_params/solar_system_atmospheres.py`.
2. Store each preset as structured data:

```python
{
    "name": "Earth",
    "planet_mass": Q_(5.972e24, "kg"),
    "planet_radius": Q_(6.371e6, "m"),
    "atmosphere_mass_fraction": 8.6e-7,
    "atomic_abundances": {...},
    "chemical_species": [...],
}
```

3. Derive atmospheric mass and column mass inside the atmosphere solver:

```text
M_atm = f_atm * M_planet
column_mass = M_atm / (4 pi R_planet^2)
```

4. Replace fixed `top_pressure` normalization with column-mass normalization for the equilibrium solver.
5. Use `z_max` only as numerical truncation and test convergence by increasing it.
6. Validate each Solar System preset against expected surface pressure order of magnitude.

## Expected Surface Pressure Checks

| Body | Expected surface pressure |
|---|---:|
| Earth | `1.0e5 Pa` |
| Mars | `~600 Pa` |
| Venus | `~9.2e6 Pa` |
| Titan | `~1.5e5 Pa` |

## Exoplanet Extension

After Solar System presets work, define class-level atmosphere mass fractions:

| Class | Suggested `f_atm` range |
|---|---:|
| airless/thin rocky | `1e-10` to `1e-7` |
| Mars-like rocky | `1e-8` to `1e-7` |
| Earth-like rocky | `1e-7` to `1e-6` |
| Venus/Titan-like thick | `1e-5` to `1e-4` |
| volatile-rich super-Earth | `1e-4` to `1e-2` |
| mini-Neptune | `1e-2` to `1e-1` |
