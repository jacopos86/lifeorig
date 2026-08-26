# ChemOrigins Local Cache

Local mirror of ChemOrigins reaction data downloaded from `https://chemorigins.bact.wisc.edu` on 2026-08-16.

## Contents

- `raw/reactions/pbr-*.json`: one ChemOrigins API JSON file per reaction.
- `raw/api_pbmdl-*.json`: ChemOrigins module API payloads.
- `raw/pbmdl-*.html`: ChemOrigins module pages.
- `raw/api_source_*.json`: ChemOrigins source API payloads where the source endpoint was available.
- `extracted/chemorigins_reactions_full.tsv`: reaction evidence table, one row per reaction-condition/source entry.
- `extracted/chemorigins_molecules_full.tsv`: molecule table extracted from reaction reactants/products.
- `extracted/chemorigins_sources_summary.tsv`: unique literature sources present in the reaction cache.
- `extracted/chemorigins_modules_summary.tsv`: ChemOrigins module summary.
- `extracted/chemorigins_all_reactions_network.txt`: parser-ready pipe network generated from all cached reactions.
- `extracted/chemorigins_environment_classification.tsv`: classification of each reaction into model subsets.
- `../PREBIOTIC/prebiotic_chemorigins_supplement_v1.txt`: ChemOrigins supplement for prebiotic core chemistry.
- `../SURFACE_POND/surface_pond_chemorigins_supplement_v1.txt`: ChemOrigins supplement for surface pond / wet-dry chemistry.
- `../VOLCANIC_ROCK_ENV/volcanic_rock_chemorigins_supplement_v1.txt`: ChemOrigins supplement for volcanic rock / mineral-surface chemistry.
- `../HYDRO_VENT/hydro_vent_chemorigins_supplement_v1.txt`: ChemOrigins supplement for hydrothermal vent chemistry.

## Counts

- Reaction JSON files: 320
- Reaction/condition evidence rows: 361
- Molecules: 240
- Unique sources: 34
- Parser-ready network reactions: 320
- Prebiotic supplement reactions: 234
- Surface pond supplement reactions: 73
- Volcanic rock supplement reactions: 73
- Hydrothermal vent supplement reactions: 75

## Notes

The generic ChemOrigins `/api/reactions` endpoint did not expose an index. The cache was built by crawling `pbr-000001` through `pbr-000750` and saving valid per-reaction API responses.

The supplement files are not mutually exclusive. A reaction can belong to both the prebiotic core and an environment supplement, for example a phosphate-compatible nucleotide precursor reaction can be both `prebiotic` and `surface_pond`.

One source API endpoint returned HTTP 500 during download:

- `10.1038/s41586-020-2330-9`

Its source metadata is still present in the per-reaction JSON files and in `extracted/chemorigins_sources_summary.tsv`.
