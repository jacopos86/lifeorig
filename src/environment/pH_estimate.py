import math
from dataclasses import dataclass
from src.common.units import Q_
from src.chemical_env.chem_env_data import ChemEnvResult
from src.utilities.logging_module import log


@dataclass
class AcidBaseSpecies:
    name: str
    kind: str                              # acid/base
    henry_ref: float                       # mol / L / atm
    eq_ref: float                          # Ka / Kb
    deltaH_henry: float | None = None      # J / Mol
    deltaH_eq: float | None = None         # J / Mol

@dataclass
class pHResult:
    pH: float
    acid_h_conc: Q_
    base_oh_conc: Q_
    net_h_conc: Q_
    species_contrib: dict[str, Q_]


ACID_BASE_SPECIES = {
    "CO2": AcidBaseSpecies(
        name="CO2",
        kind="acid",
        henry_ref=3.3e-2,
        eq_ref=4.45e-7,
        deltaH_henry=None,
        deltaH_eq=None
    ),
    "NH3": AcidBaseSpecies(
        name="NH3",
        kind="base",
        henry_ref=5.9e1,
        eq_ref=1.8e-5,
        deltaH_henry=None,
        deltaH_eq=None
    ),
    "H2S": AcidBaseSpecies(
        name="H2S",
        kind="acid",
        henry_ref=1.0e-1,
        eq_ref=9.1e-8,
        deltaH_henry=None,
        deltaH_eq=None
    ),
}


def estimate_water_pH_from_chemistry(
        chemical_env: ChemEnvResult,
        species_data: dict[str, AcidBaseSpecies] = ACID_BASE_SPECIES,
    ) -> pHResult:
    if chemical_env is None:
        log.error("chemical_env is required to estimate pH")

    acid_h_conc = Q_(0.0, "mol / liter")
    base_oh_conc = Q_(0.0, "mol / liter")
    species_contrib = {}

    for species, data in species_data.items():
        partial_pressure = chemical_env.get_partial_pressure(species, default=None)
        if partial_pressure is None:
            continue
        partial_pressure_atm = partial_pressure.to("atm").magnitude
        if partial_pressure_atm <= 0.0:
            continue

        dissolved = Q_(data.henry_ref * partial_pressure_atm, "mol / liter")
        contribution = Q_(
            math.sqrt(max(data.eq_ref * dissolved.to("mol / liter").magnitude, 0.0)),
            "mol / liter"
        )
        species_contrib[species] = contribution
        if data.kind == "acid":
            acid_h_conc += contribution
        elif data.kind == "base":
            base_oh_conc += contribution
        else:
            log.error(f"Unknown acid/base species kind: {data.kind}")

    water_neutral = Q_(1.0e-7, "mol / liter")
    net_h_conc = acid_h_conc - base_oh_conc + water_neutral
    if net_h_conc.magnitude <= 0.0:
        net_h_conc = water_neutral

    pH = -math.log10(net_h_conc.to("mol / liter").magnitude)
    return pHResult(
        pH=pH,
        acid_h_conc=acid_h_conc,
        base_oh_conc=base_oh_conc,
        net_h_conc=net_h_conc,
        species_contrib=species_contrib,
    )
