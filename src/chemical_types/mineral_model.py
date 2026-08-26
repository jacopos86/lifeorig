from dataclasses import dataclass, field

MINERAL_PHASE_TAGS = ("(aq)", "(aer)", "(s)", "(sol)")

KNOWN_MINERAL_NAMES = {
    "apatite",
    "basalt_glass",
    "brucite",
    "carbonate_surface",
    "clay",
    "FeS",
    "FeS2",
    "FeO_silicate",
    "FeOOH",
    "Fe_oxide",
    "magnetite",
    "mineral_surface",
    "NiS",
    "olivine",
    "phosphate_surface",
    "pyrite",
    "pyroxene",
    "serpentine_like_surface",
    "silica",
    "silica_surface",
}

@dataclass
class MineralSpecies:
    ID: int
    name: str
    aliases: set[str]
    phase: str
    catalytic_activity: dict[str, float] = field(default_factory=dict)

@dataclass
class MineralTemplate:
    ID: int
    name: str
    aliases: set[str]
    phase: str
    template_type: str
    matching_minerals: list[MineralSpecies]

def mineral_base_name(species_name):
    for phase_tag in MINERAL_PHASE_TAGS:
        if species_name.endswith(phase_tag):
            return species_name[: -len(phase_tag)]
    return species_name

def is_mineral(species_name):
    return mineral_base_name(species_name) in KNOWN_MINERAL_NAMES
