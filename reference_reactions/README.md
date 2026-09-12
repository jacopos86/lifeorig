# LIFEORIG Reaction-Network Catalogue

This directory contains reaction data for different physical environments.
Files from different environments must not be merged automatically: atmospheric
gas kinetics, aqueous chemistry, mineral-surface chemistry, and polymerization
use different rate laws and state variables.

## File categories

- **Raw source**: an original downloaded model/database file. Preserve it; parse
  or convert it before use.
- **Full kinetic network**: suitable as the primary mechanism for its stated
  environment when its rate laws are implemented.
- **Selected subset**: useful for tests or a restricted module, but incomplete.
- **Combined network**: assembled from several sources and deduplicated.
- **Candidate/hypothesis network**: useful for model development, not a
  validated kinetic mechanism.
- **Supplement/audit file**: provenance or additional records, not normally a
  standalone simulation network.

## Atmospheric networks

### Early Earth

| File | Category | Contents | Recommended use |
| --- | --- | --- | --- |
| `ATMOSPHERIC/VPL_ATMOS/PHOTOCHEM/INPUTFILES/TEMPLATES/Archean+haze/reactions.rx` | Full kinetic network, raw VPL format | 392 Archean atmospheric gas-phase and photolysis reactions with rate parameters | Primary network for `NonEquilLayeredAtmosphSolver` |
| `EARTH_LIKE/earth_like_atmospheric_photochemistry.txt` | Selected subset | 20 illustrative photolysis/radical reactions with mostly symbolic rates | Small parser/debugging tests only |
| `EARTH_LIKE/earth_like_hadean_photolysis_selected.txt` | Selected subset | Photolysis branches extracted from the VPL Archean mechanism | Photolysis-only tests; do not use instead of the full mechanism |

The full VPL file is documented beside the source in
`ATMOSPHERIC/VPL_ATMOS/PHOTOCHEM/INPUTFILES/TEMPLATES/Archean+haze/README.md`.

### Titan

| File | Category | Contents | Recommended use |
| --- | --- | --- | --- |
| `TITAN/titan_gas_phase_hebrard2013.txt` | Full kinetic network, converted | Complete downloaded KIDA Hébrard 2013 gas-phase rows in LIFEORIG pipe format | Primary Titan gas-phase chemistry |
| `TITAN/titan_photolysis_selected.txt` | Selected/combined subset | Titan photolysis branches assembled from candidate and VPL sources | Add photolysis to the KIDA gas network |
| `TITAN/titan_hebrard2013_bimolecular_selected.txt` | Selected subset | Small selection from KIDA bimolecular reactions | Tests and comparison only |
| `TITAN/titan_vpl_photolysis_selected.txt` | Selected subset | VPL Titan photolysis branches | Provenance/comparison; largely incorporated into `titan_photolysis_selected.txt` |
| `TITAN/titan_reaction_network_hebrard2013_combined.txt` | Combined network | KIDA gas reactions plus nonduplicate candidate Titan reactions | Broad single-file Titan candidate mechanism |
| `TITAN/titan_reaction_network_v0.txt` | Candidate network | Literature-guided atmospheric scaffold | Historical/development file; rates require validation |
| `TITAN/titan_reactions_by_environment_v0.txt` | Candidate multi-environment network | Gas, lake, solid, aerosol and transient aqueous reactions | Environment-classification design; not one homogeneous kinetic solver |
| `TITAN/titan_liquid_surface_polymerization.txt` | Selected environment module | Liquid, aerosol and surface polymerization reactions | Titan surface/liquid solver only |
| `TITAN/titan_polymerization_catalysis_module_v0.txt` | Candidate module | Original polymerization/catalysis scaffold | Historical source for the converted liquid/surface file |

Raw KIDA source data are preserved under `ATMOSPHERIC/KIDA/Hebrard2013/`.

## Earth-like surface and prebiotic networks

| File | Category | Contents | Recommended use |
| --- | --- | --- | --- |
| `EARTH_LIKE/earth_like_prebiotic_evolution_network.txt` | Selected environment module | HCN-derived prebiotic evolution reactions without hydrothermal feedstock blocks | Earth-like surface/prebiotic evolution, not atmospheric kinetics |
| `EARTH_LIKE/earth_like_volcanic_rock_combined_network_v1.txt` | Combined network | Deduplicated photolysis, volcanic-rock and prebiotic reactions | Coupled volcanic-rock environment after phase filtering |
| `EARTH_LIKE/earth_like_volcanic_rock_combined_duplicates_v1.tsv` | Audit file | Duplicate-removal decisions for the combined network | Provenance only; never load as a reaction network |
| `VOLCANIC_ROCK_ENV/volcanic_rock_reaction_network_v1.txt` | Candidate environment network | Gas-rock, mineral-surface, porous-film and wet/dry volcanic chemistry | Volcanic-rock environment solver |
| `VOLCANIC_ROCK_ENV/volcanic_rock_chemorigins_supplement_v1.txt` | Supplement | ChemOrigins experimental records classified as volcanic/mineral chemistry | Add evidence-backed reactions selectively |
| `HYDRO_VENT/hydro_vent_reaction_network_v1.txt` | Candidate environment network | High-pressure aqueous, serpentinization and vent-ocean chemistry | Hydrothermal-vent solver |
| `HYDRO_VENT/hydro_vent_chemorigins_supplement_v1.txt` | Supplement | ChemOrigins hydrothermal experimental records | Add selectively to the vent network |
| `SURFACE_POND/surface_pond_reaction_network_v1.txt` | Candidate environment network | UV, wet/dry, eutectic and evaporative pond chemistry | Surface-pond solver |
| `SURFACE_POND/surface_pond_chemorigins_supplement_v1.txt` | Supplement | ChemOrigins pond/prebiotic experimental records | Add selectively to the pond network |

## General prebiotic networks

| File | Category | Contents | Recommended use |
| --- | --- | --- | --- |
| `PREBIOTIC/prebiotic_cyanosulfidic_literature_network_v1.txt` | Literature-backed candidate network | HCN/cyanosulfidic pathways with explicit literature references | Preferred compact prebiotic chemistry module |
| `PREBIOTIC/prebiotic_polymer_autocatalytic_network_v2.txt` | Hypothesis network | Rule-based polymer growth and autocatalytic closure | Autocatalysis experiments, not validated chemistry |
| `PREBIOTIC/hcn_thioester_autocatalytic_network_v0.txt` | Historical candidate network | Early HCN/thioester scaffold | Retained for provenance; prefer later networks |
| `PREBIOTIC/hcn_thioester_autocatalytic_network_v0.csv` | Historical parameter table | Original scaffold with unresolved rate placeholders | Provenance only |
| `PREBIOTIC/hcn_thioester_autocatalytic_network_v1.csv` | Parameterized candidate table | Development rates and catalysis probabilities | Algorithm tests, not literature kinetics |
| `PREBIOTIC/prebiotic_chemorigins_supplement_v1.txt` | Supplement | General ChemOrigins experimental records | Select by environment before use |
| `PREBIOTIC/chemorigins_extracted_records_v1.tsv` | Audit/cache file | Manually extracted ChemOrigins provenance records | Provenance only |

`ChemOrigins/` contains the larger raw and extracted source dataset. See its own
`README.md`; do not treat the full export as one physically consistent network.

## Current solver selection

The non-equilibrium early-Earth atmosphere currently selects only:

```text
ATMOSPHERIC/VPL_ATMOS/PHOTOCHEM/INPUTFILES/TEMPLATES/Archean+haze/reactions.rx
```

Surface and volcanic-rock networks remain separate. Coupling between them
should occur through boundary fluxes such as deposition, rainout, escape,
outgassing, and atmosphere-ocean exchange.
