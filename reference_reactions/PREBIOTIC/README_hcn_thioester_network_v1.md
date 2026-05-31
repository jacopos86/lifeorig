# HCN/Thioester Autocatalytic Network v1

This folder contains a tentative prebiotic reaction network intended for
autocatalysis threshold tests.

Use this file for computation:

```text
hcn_thioester_autocatalytic_network_v1.csv
```

The older `.txt` file is human-readable notes. The `.csv` file is the
machine-facing version.

## CSV Columns

- `reaction_id`: stable reaction identifier.
- `section`: chemical subsystem, for example `HCN_core`,
  `strecker_amino_acids`, `thioester_energy`, or `autocatalytic_closure`.
- `reaction_type`: qualitative class of reaction.
- `reactants`: explicit reactant species. Species are separated by `;`.
- `products`: explicit product species. Species are separated by `;`.
- `catalysts`: species or surfaces that catalyze the reaction. If a catalyst is
  also present in reactants/products, it is still listed here for rate-law use.
- `conditions`: environmental condition set, e.g. `K02`, `K06`, `K10`.
- `rate_law`: currently `mass_action`.
- `k_uncatalyzed`: baseline rate coefficient placeholder.
- `k_catalyzed`: catalyzed rate coefficient placeholder.
- `catalysis_probability`: probability that the listed catalyst is active in a
  stochastic/protocell-network simulation.
- `autocatalytic_role`: why the reaction is included.
- `source_hint`: literature/mechanism family hint.

## Current Rate Convention

For a reaction with reactants `A; B` and catalyst `C`, use:

```text
rate = k_uncatalyzed [A][B]
     + catalysis_probability * k_catalyzed [A][B][C]
```

For photolysis reactions with `hnu`, replace the mass-action expression with a
photolysis rate later:

```text
rate = J_species [species]
```

The numeric values in `k_uncatalyzed`, `k_catalyzed`, and
`catalysis_probability` are placeholders for threshold experiments. They are
not literature rate constants.

## Important Caveat

This is not a validated biochemical network. It is a structured hypothesis
network designed to test whether feedstock fluxes plus catalytic probability
can cross an autocatalytic threshold.
