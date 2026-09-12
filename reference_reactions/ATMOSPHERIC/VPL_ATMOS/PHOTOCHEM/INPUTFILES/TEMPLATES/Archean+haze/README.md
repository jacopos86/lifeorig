# VPL Archean Atmospheric Reaction Network

## Reaction file

`reactions.rx` is the full atmospheric chemical mechanism used by the VPL
Atmos Archean-Earth haze template.

The mechanism contains 392 reactions involving approximately 76 atmospheric
species:

- 257 bimolecular reactions (`2BODY`)
- 15 termolecular and pressure-dependent reactions (`3BODY`)
- 60 reactions with specialized rate expressions (`WEIRD`)
- 60 photolysis reactions (`PHOTO`, `PHOTP`, or `PHOTX`)

This is an atmospheric gas-phase and photochemical network. It is not an
aqueous, mineral-surface, or volcanic-rock reaction network.

## LIFEORIG usage

The file is selected in an input configuration relative to
`reference_reactions`:

```json
{
    "chemical_network": {
        "type": "atmospheric",
        "reaction_files": [
            "ATMOSPHERIC/VPL_ATMOS/PHOTOCHEM/INPUTFILES/TEMPLATES/Archean+haze/reactions.rx"
        ]
    }
}
```

`VPLAtmosReactionParser` parses the native VPL format. The reaction database
driver then registers the parsed reactions and their source-file hash in
MySQL. `NonEquilLayeredAtmosphSolver` reads the selected network back from the
database.

`HV` in the source file is normalized to `hnu`. Photons and the third-body
symbol `M` remain reaction controls and must not be integrated as chemical
species.

## Rate format

The source preserves the original VPL reaction classes and rate data. These
include Arrhenius-like coefficients, falloff parameters, explicit
temperature/density expressions, and photolysis branches. LIFEORIG currently
stores these expressions verbatim; evaluation by the non-equilibrium kinetic
solver is a separate implementation step.

## References

- Arney, G. et al. (2016), *The Pale Orange Dot: The Spectrum and
  Habitability of Hazy Archean Earth*, Astrobiology 16, 873-899.
  DOI: 10.1089/ast.2015.1422.
- Zerkle, A. L. et al. (2012), *A bistable organic-rich atmosphere on the
  Neoarchaean Earth*, Nature Geoscience 5, 359-363.
  DOI: 10.1038/ngeo1425.
- Virtual Planetary Laboratory Atmos model:
  https://github.com/VirtualPlanetaryLaboratory/atmos

Individual reaction comments cite the underlying kinetic evaluations where
available, including JPL and NIST sources.
