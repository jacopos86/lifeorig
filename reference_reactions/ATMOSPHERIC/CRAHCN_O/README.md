# CRAHCN atmospheric networks

This directory contains two separately importable gas-phase kinetic sources. They are intentionally not merged with the VPL Atmos file. LIFEORIG loads the files in the configured order and removes duplicate net equations in memory, retaining the first occurrence.

## `crahcn_o_pearce2020.txt`

This is the complete 104-reaction CRAHCN network followed by its published oxygen extension. Multi-step mechanisms in the supporting tables are represented by their net reactants and final products. Electronic-state names are mapped to the VPL convention (`O1D`, `O`, `CH21`, `CH23`, `N2D`, and `N`).

Sources:

- Pearce et al. (2020), *HCN Production in Titan's Atmosphere: Coupling Quantum Chemistry and Disequilibrium Atmospheric Modeling*, DOI: 10.3847/1538-4357/abae5b.
- Pearce et al. (2020), *An Experimental and Theoretical Investigation of HCN Production in the Hadean Earth Atmosphere*, DOI: 10.1021/acs.jpca.0c06804.

Rate laws are copied from the published Lindemann and modified-Arrhenius tables. CRAHCN and CRAHCN-O are reduced networks designed primarily for HCN and formaldehyde chemistry over approximately 50–400 K; they are not universal combustion networks.

## `hadean_nh3_no_extension_pearce2022.txt`

This file contains the 53 added two-body reactions from Table S1 of Pearce, He, and Hörst (2022). It expands hydrocarbon chemistry and explicitly adds NH3 and NO chemistry.

Source: Pearce, He, and Hörst (2022), *An Experimental and Theoretical Investigation of HCN Production in the Hadean Earth Atmosphere*, DOI: 10.1021/acsearthspacechem.2c00138.

The paper also removed 12 older CRAHCN-O reactions whose products lacked sinks. The supporting table identifies the 53 additions but does not enumerate all 12 removals. Consequently, this file is kept as an additive source, and those removals must not be guessed silently.

## Reproduction

The files were transcribed from the authors' arXiv LaTeX supporting tables using `scripts/import_crahcn_literature.py`. The importer preserves reaction IDs, rate parameters, paper/table provenance, and keeps the two published sources separate.
