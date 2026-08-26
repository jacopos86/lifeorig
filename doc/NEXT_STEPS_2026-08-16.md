# Work Memo - 2026-08-16

Goal: recover the state of the Titan / Earth-like reaction-network work and define the next steps.

## Current Direction

The code is now being pushed toward a workflow where each planet/environment scenario uses one deduplicated reaction file:

- Titan scenario: Titan photochemistry + gas reactions + liquid/surface/polymerization modules.
- Earth-like volcanic rock scenario: Hadean photolysis + volcanic rock + prebiotic chemistry + ChemOrigins supplements.
- Later Earth-like scenarios: same structure for surface pond and hydrothermal vent.

The molecular classifier is becoming the central object builder. It must turn reaction species names into unique entities with aliases, phases, states, conformations, templates, minerals, and generated polymer members.

## Earth-Like Current State

The active Earth-like script is:

```text
TEST_SCRIPTS/run_earth_like.sh
```

It now reads a single combined deduplicated file:

```text
reference_reactions/EARTH_LIKE/earth_like_volcanic_rock_combined_network_v1.txt
```

This file includes:

- Hadean/Archean photolysis extracted from VPL `Archean+haze/reactions.rx`.
- Volcanic rock reaction network.
- Literature-backed prebiotic cyanosulfidic network.
- ChemOrigins prebiotic supplement.
- ChemOrigins volcanic-rock/mineral supplement.

Duplicate audit:

```text
reference_reactions/EARTH_LIKE/earth_like_volcanic_rock_combined_duplicates_v1.tsv
```

Current parsed count in `earth_like_volcanic_rock_combined_network_v1.txt`:

```text
516 reactions
501 species
```

Important: the script uses `environment_source="explicit"`, so it should not run the atmospheric/planetary solver.

Run command:

```bash
LIFEORIG_REACTION_DB_PASSWORD='...' TEST_SCRIPTS/run_earth_like.sh
```

## ChemOrigins Cache

Full local ChemOrigins cache was created:

```text
reference_reactions/ChemOrigins/
```

Key files:

```text
reference_reactions/ChemOrigins/raw/reactions/
reference_reactions/ChemOrigins/extracted/chemorigins_all_reactions_network.txt
reference_reactions/ChemOrigins/extracted/chemorigins_reactions_full.tsv
reference_reactions/ChemOrigins/extracted/chemorigins_molecules_full.tsv
reference_reactions/ChemOrigins/extracted/chemorigins_sources_summary.tsv
reference_reactions/ChemOrigins/extracted/chemorigins_environment_classification.tsv
```

Counts:

```text
320 reaction JSON files
361 reaction/source/condition evidence rows
240 molecules
34 literature sources
```

Supplement files generated:

```text
reference_reactions/PREBIOTIC/prebiotic_chemorigins_supplement_v1.txt
reference_reactions/SURFACE_POND/surface_pond_chemorigins_supplement_v1.txt
reference_reactions/VOLCANIC_ROCK_ENV/volcanic_rock_chemorigins_supplement_v1.txt
reference_reactions/HYDRO_VENT/hydro_vent_chemorigins_supplement_v1.txt
```

These supplements are not mutually exclusive. A reaction can be both prebiotic and environment-relevant.

## Titan Current State

The active Titan script is:

```text
TEST_SCRIPTS/run_titan.sh
```

It currently reads:

```text
reference_reactions/TITAN/titan_reaction_network_v0.txt
```

Titan has a stronger explicit polymerization structure than Earth-like:

```text
P_HCN_n
P_C2H2_n
P_CN_n
P_CxHy_n
PAH_seed
PAH_large
N_PAH
tholin / residue templates
```

Important Titan files:

```text
reference_reactions/TITAN/titan_liquid_surface_polymerization.txt
reference_reactions/TITAN/titan_polymerization_catalysis_module_v0.txt
reference_reactions/TITAN/titan_vpl_photolysis_selected.txt
reference_reactions/TITAN/titan_hebrard2013_bimolecular_selected.txt
```

Need to finish validating whether `titan_reaction_network_v0.txt` is the correct final deduplicated file, or whether Titan should also be rebuilt into one final combined deduplicated network from its components.

Run command:

```bash
LIFEORIG_REACTION_DB_PASSWORD='...' TEST_SCRIPTS/run_titan.sh
```

## Molecular Classifier State

Main files:

```text
src/chemical_types/define_molecule_set.py
src/chemical_types/molecule_model.py
src/chemical_types/mineral_model.py
```

Current important ideas:

- `ReferenceMolecule` should represent a unique molecular entity with ID, sequence, aliases, phase, state, conformation.
- Aliases must include all reaction names that refer to the same entity.
- Templates do not own independent matrix IDs; their matching molecules do.
- Mineral species have catalytic activity.
- `mineral*` is a mineral template.
- Atoms are handled separately from molecules.
- Phase should be per entity; same sequence in different phase can be separate if needed for network entries.
- Excited states are separate entities when needed.
- Conformation matters for `cC3H`, `lC3H`, `tC3H2`, and later for polymers.

Known classifier questions still open:

- How to define generated peptide and oligonucleotide templates cleanly.
- How to map ChemOrigins long names and `pbm-*` placeholders into stable molecule entities.
- How to expand generic rules such as `peptide_n + amino_acid -> peptide_n+1`.
- How to assign generated polymers to templates and matrix IDs.
- How to keep aliases complete without merging true isomers incorrectly.

## Polymerization And Catalysis

Titan already has explicit polymer family chemistry.

Earth-like volcanic rock has generic polymerization / oligomerization proxies:

```text
amino_acid + amino_acid -> peptide2 + H2O
peptide_n + amino_acid -> peptide_n+1 + H2O
nucleotide + nucleotide -> dinucleotide + H2O
oligonucleotide_n + nucleotide -> oligonucleotide_n+1 + H2O
template_oligo + nucleotide -> extended_template
oligo_n + monomer -> oligo_n+1
```

Earth-like also has hydrolysis/reverse/loss reactions:

```text
peptide_n + H2O -> peptide_n-1 + amino_acid
oligonucleotide_n + H2O -> oligonucleotide_n-1 + nucleotide
```

Photolysis exists mainly for atmospheric/small molecules. Polymer photolysis fragmentation is not yet explicit:

```text
peptide_n + hnu -> peptide_fragments
oligonucleotide_n + hnu -> oligo_fragments
HCN_polymer_n + hnu -> nitrile_fragments
```

Future catalysis model:

```text
rate = base_rate * environment_factor * catalyst_factor
```

`catalyst_factor` should come from generated polymer properties:

- sequence
- length
- conformation/folding
- phase/environment
- surface binding
- hydrophobic/polar pattern
- reactant class specificity

## Planetary Solver Issue

`src/__main__.py` used to run `planetary_solver_driver(p)` unconditionally for every `set_initial_state`.

This was changed so it only runs when:

```python
p.environment_config.uses_planetary_solver
```

For explicit Earth-like volcanic rock runs:

```json
"environment_source": "explicit"
```

so the code logs:

```text
SKIPPING PLANETARY SOLVER: environment_source is explicit
```

This needs cleanup later. The current naming is functional but not final. We should eventually separate:

- planetary metadata / preset planet object
- atmospheric equilibrium solver
- local environment builder
- reaction-network import

## Immediate Next Steps

1. Run Earth-like volcanic rock with DB password and confirm it reaches molecular classifier.
2. Fix classifier failures from the 501-species Earth-like combined network.
3. Finish Titan classifier validation with the Titan reaction set.
4. Build deduplicated combined files for:

```text
EARTH_LIKE/earth_like_surface_pond_combined_network_v1.txt
EARTH_LIKE/earth_like_hydro_vent_combined_network_v1.txt
```

5. Decide whether Titan needs a regenerated deduplicated combined file like Earth-like.
6. Clean up `uses_planetary_solver` naming and control flow.
7. Return to the planetary/atmospheric solver after the reaction-network and molecular classifier are stable.

## Do Not Forget

- Do not include component reaction files together with a combined file in a run script.
- One final deduplicated file per environment scenario.
- Keep ChemOrigins raw data intact.
- Use duplicate audit files to decide if a reaction should really be merged or kept separate.
- Hydrothermal reactions should only enter hydro vent scenarios unless explicitly modeling transport between habitats.
