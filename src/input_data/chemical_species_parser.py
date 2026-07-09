import re
from dataclasses import dataclass


@dataclass(frozen=True)
class ChemicalSpecies:
    name: str


@dataclass(frozen=True)
class Molecule(ChemicalSpecies):
    phase: str | None = None


@dataclass(frozen=True)
class IonicSpecies(ChemicalSpecies):
    charge: int
    phase: str | None = None


@dataclass(frozen=True)
class Mineral(ChemicalSpecies):
    phase: str = "solid"


@dataclass(frozen=True)
class UncategorizedMolecularSpecies(ChemicalSpecies):
    reason: str = "unclassified"


@dataclass(frozen=True)
class MineralSurface(ChemicalSpecies):
    surface_type: str | None = None


KNOWN_MINERALS = {
    "apatite",
    "FeOOH",
    "FeS",
    "Mg_phosphate_solid",
    "Fe_phosphate_solid",
    "Mg_fatty_acid_solid",
    "Ca_fatty_acid_solid",
}

MOLECULE_NAMES = {
    "acetaldehyde",
    "acetamide",
    "acetate",
    "acetic_acid",
    "acetyl_phosphate",
    "acetyl_SR",
    "acyl_phosphate_ester",
    "alanine",
    "alanine_nitrile",
    "alpha_ketoglutarate",
    "amidine_pool",
    "amino_acid",
    "amino_acid_pool",
    "aminoacetonitrile",
    "aminoacyl_phosphate",
    "aminoacyl_thioester",
    "aminoamide_pool",
    "aminonitrile_pool",
    "amino_oxazole_derivative",
    "aspartate",
    "carbamate",
    "CH3COOH",
    "CH3COSR",
    "CH3OH",
    "CH3SH",
    "CH4",
    "CO",
    "CO2",
    "COS",
    "cyanamide",
    "fatty_acid_Cn",
    "formaldehyde",
    "formaldimine",
    "formamide",
    "formate",
    "formic_acid",
    "fumarate",
    "glutamate",
    "glycerol_like_alcohol",
    "glycine",
    "glycolaldehyde",
    "glycolate",
    "glycolonitrile",
    "glyoxylate",
    "H2",
    "H2O",
    "H2S",
    "HCN",
    "HCHO",
    "HCONH2",
    "hetero_oligo_n",
    "hydroxymethylene_carbene",
    "imine_pool",
    "imino_diacetate",
    "isocyanic_acid",
    "lactate",
    "malate",
    "metallopeptide_complex",
    "micelle",
    "monomer_pool",
    "N_acetyl_amino_acid",
    "N_acyl_amino_acid",
    "N2",
    "NH2CO-SH",
    "NH3",
    "organic_degradation_pool",
    "organic_pool",
    "oxidized_products",
    "peptide2",
    "peptide_m",
    "peptide_m+n",
    "peptide_n",
    "peptide_n+1",
    "peptide_n-1",
    "peptide_vesicle",
    "protonated_fatty_acid_Cn",
    "pyrophosphate",
    "pyruvate",
    "redox_products",
    "refractory_organic_surface",
    "RSH",
    "smaller_hetero_oligo",
    "succinate",
    "thioacetamide",
    "thioacetate",
    "thiocarbamate",
    "thioformamide",
    "thiophosphate",
    "urea",
    "vesicle",
    "vesicle_fragments",
}


def parse_species_into_chemicals(species_names):
    return [classify_species(name) for name in species_names]


def classify_species(name):
    phase = parse_phase(name)
    if is_mineral_surface(name):
        return MineralSurface(name=name, surface_type=parse_surface_type(name))
    if is_mineral(name, phase):
        return Mineral(name=name)
    charge = parse_charge(name)
    if charge is not None:
        return IonicSpecies(name=name, charge=charge, phase=phase)
    if phase == "aqueous" or name in MOLECULE_NAMES:
        return Molecule(name=name, phase=phase)
    return UncategorizedMolecularSpecies(name=name)


def parse_phase(name):
    if name.endswith("(aq)"):
        return "aqueous"
    if name.endswith("(s)") or name.endswith("(silicate)"):
        return "solid"
    return None


def is_mineral_surface(name):
    lowered = name.lower()
    return (
        "surface" in lowered
        or lowered.endswith("_patch")
        or lowered.endswith("_interface")
    )


def parse_surface_type(name):
    if "surface" in name:
        return name.split("surface", 1)[0].rstrip("_") or "generic"
    if name.endswith("_patch"):
        return "patch"
    if name.endswith("_interface"):
        return "interface"
    return None


def is_mineral(name, phase):
    return phase == "solid" or name in KNOWN_MINERALS or name.endswith("_solid")


def parse_charge(name):
    if "_" in name:
        return None
    match = re.search(r"([+-]+)$", name)
    if match is None:
        return None
    signs = match.group(1)
    if signs.startswith("+"):
        return len(signs)
    return -len(signs)
