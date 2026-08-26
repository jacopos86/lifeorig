import re
from itertools import product
from dataclasses import dataclass
from periodictable import elements
from src.chemical_types.molecule_model import MultiPolymer
from src.chemical_types.molecule_model import ReferenceMolecule, MolecularTemplate, AtomicSpecies
from src.chemical_types.mineral_model import MineralSpecies, MineralTemplate, is_mineral
from src.network_generation.reaction_mysql_db import open_reaction_database, get_species_from_reaction_database
from src.utilities.logging_module import log

NON_MOLECULAR_SPECIES = {
    "M",
    "hnu",
    "hv",
    "photon",
}

NAMED_MOLECULE_SEQUENCES = {
    "cyanoacetamide": "NCCH2CONH2",
    "formamide_HCN_adduct": "HCONH2.HCN",
}

TEMPLATE_POLYMER_FAMILIES = {
    "C_rich_polymer": ("C2H2",),
    "N_rich_polymer": ("HCN", "CN"),
    "N_rich_soluble_polymer": ("HCN", "CN"),
    "P_C2H2_n": ("C2H2",),
    "P_C2H2_nplus1": ("C2H2",),
    "P_CN_n": ("CN",),
    "P_CN_nplus1": ("CN",),
    "P_HCN_n": ("HCN",),
    "P_HCN_nplus1": ("HCN",),
}

PAH_TEMPLATE_NAMES = {
    "PAH_seed",
    "PAH_large",
    "N_PAH",
}

PAH_TEMPLATE_REPRESENTATIVES = {
    "PAH_seed": (
        ("C10H8", "naphthalene"),
    ),
    "PAH_large": (
        ("C14H10", "anthracene"),
        ("C14H10", "phenanthrene"),
    ),
    "N_PAH": (
        ("C15H9N", "cyanoanthracene"),
        ("C15H9N", "cyanophenanthrene"),
    ),
}

#
#   build the molecules set
#

def build_molecular_string_model_set(metabolites_data):
    """
    Build the global set of unique Molecule objects present in the simulation.
    Returns:
        molecules : list[Molecule]
        molecule_map : dict[sequence -> Molecule]
    """
    max_polymer_length = metabolites_data.get("pol_strng_maxsize")
    mol_type = metabolites_data.get("type")
    molecules = []
    molecule_map = {}
    if mol_type == "binary":
        # build all binary polymers up to max length
        for L in range(1, max_polymer_length + 1):
            for i in range(2 ** L):
                seq = format(i, f'0{L}b')  # binary string
                #mol = Molecule(seq)
                molecules.append(mol)
                molecule_map[seq] = mol
    elif mol_type == 'multi':
        alphabet = list(MultiPolymer._MONOMER_TYPES)
        for L in range(1, max_polymer_length + 1):
            for seq in product(alphabet, repeat=L):
                #mol = Molecule(list(seq))
                molecules.append(mol)
                molecule_map[tuple(seq)] = mol
    else:
        log.error("Unknown molecule type")
    log.info("\n")
    for i, mol in enumerate(molecules):
        log.info(f"\t molecule {i+1}: {mol.show_sequence()} -- length: {len(mol)}")
    log.info("\t " + p.sep)
    log.info("\n")
    return molecules, molecule_map

#
#   build complete molecular species set
#

def build_molecular_species_set(input_params, reaction_source_files=None):
    species = set()
    # collect chemical species
    species.update(_species_from_reaction_database(reaction_source_files))
    species.update(_species_from_planetary_chemistry(input_params.planetary_data, input_params.env_model))
    species.update(_species_from_local_environment(input_params.local_env_data))
    max_polymer_length = (input_params.metabolites_params or {}).get("pol_strng_maxsize")
    return MolecularSpeciesSet.from_species_names(species, max_polymer_length=max_polymer_length)

#
#   get data from reaction file
#

def _species_from_reaction_database(reaction_source_files=None):
    db = open_reaction_database()
    try:
        return set(get_species_from_reaction_database(db, reaction_source_files))
    finally:
        db.close()

#
#   get species from planet chemistry
#

def _species_from_planetary_chemistry(planet, env_model):
    if planet is None:
        return set()
    if env_model == "hydro_vent":
        hydro = planet.hydro or {}
        return set((hydro.get("ocean_composition") or {}).keys())
    chemistry = planet.chemical_stationary_config
    if chemistry is not None:
        return set(chemistry.atmosph_composition(default={}).keys())
    chemistry_input = planet.chemistry or {}
    return set(chemistry_input.get("chemical_species") or [])

#
#   get solvent species
#

def _species_from_local_environment(local_env_data):
    if not local_env_data:
        return set()
    solvent_data = local_env_data.get("solvent_data")
    if solvent_data is None:
        return set()
    return set(solvent_data.composition.keys())

#
#     Molecular species set construction
#

@dataclass
class MolecularSpeciesSet:
    species: list

    PHASE_TAGS = {
        "(aer)": "aerosol",
        "(aq)": "aqueous",
        "(s)": "solid",
        "(sol)": "solution",
    }

    @staticmethod
    def _atomic_symbols():
        return {
            element.symbol
            for element in elements
            if getattr(element, "symbol", None)
        }

    @classmethod
    def _split_phase(cls, species_name):
        for phase_tag, phase in cls.PHASE_TAGS.items():
            if species_name.endswith(phase_tag):
                return species_name[: -len(phase_tag)], phase
        return species_name, "gas"

    @staticmethod
    def _split_state(species_name):
        species_name = re.sub(r"(?<!\d)1(?!\d)", "", species_name)
        if species_name.endswith("X"):
            return _canonical_formula(species_name[:-1]), "excited"
        if _is_numbered_polymer_species(species_name):
            return _canonical_polymer_formula(species_name), "ground"
        return _canonical_formula(species_name), "ground"

    @staticmethod
    def _split_conformation(species_name):
        conformation_prefixes = {
            "c": "cyclic",
            "l": "linear",
            "t": "trans",
        }
        if len(species_name) > 1 and species_name[0] in conformation_prefixes and species_name[1].isupper():
            return species_name[1:], conformation_prefixes[species_name[0]]
        return species_name, "linear"

    @staticmethod
    def _split_atomic_state(species_name):
        match = re.fullmatch(r"([A-Z][a-z]?)(?:\(([^)]+)\))?", species_name)
        if match is None:
            return species_name, "ground"
        symbol, state = match.groups()
        if state == "4S":
            state = "ground"
        return symbol, state or "ground"

    @classmethod
    def from_species_names(cls, species_names, max_polymer_length=None):
        atomic_aliases = {}
        molecule_aliases = {}
        template_aliases = {}
        mineral_aliases = {}
        mineral_template_aliases = {}
        atomic_symbols = cls._atomic_symbols()
        for species_name in sorted(species_names):
            print(species_name)
            if species_name in NON_MOLECULAR_SPECIES:
                continue
            species_base_name, phase = cls._split_phase(species_name)
            atomic_symbol, atomic_state = cls._split_atomic_state(species_base_name)
            if atomic_symbol in atomic_symbols and species_base_name == atomic_symbol:
                atomic_key = (atomic_symbol, phase, atomic_state)
                atomic_aliases.setdefault(atomic_key, set()).add(species_name)
                continue
            if atomic_symbol in atomic_symbols and species_base_name != atomic_symbol:
                atomic_key = (atomic_symbol, phase, atomic_state)
                atomic_aliases.setdefault(atomic_key, set()).add(species_name)
                continue
            if _is_mineral_template(species_base_name):
                mineral_template_key = (species_base_name, phase)
                mineral_template_aliases.setdefault(mineral_template_key, set()).add(species_name)
            elif _is_mineral_species(species_name):
                mineral_key = (species_base_name, phase)
                mineral_aliases.setdefault(mineral_key, set()).add(species_name)
            elif _is_molecular_template(species_base_name):
                template_key = (species_base_name, phase, "ground")
                template_aliases.setdefault(template_key, set()).add(species_name)
            else:
                species_base_name, state = cls._split_state(species_base_name)
                species_base_name, conformation = cls._split_conformation(species_base_name)
                molecule_key = (species_base_name, phase, state, conformation)
                molecule_aliases.setdefault(molecule_key, set()).add(species_name)

        max_polymer_length = _resolve_max_polymer_length(max_polymer_length, molecule_aliases)
        _add_template_polymer_molecules(template_aliases, molecule_aliases, max_polymer_length)

        species = []
        next_id = 1
        for atomic_symbol, phase, state in sorted(atomic_aliases):
            species.append(
                AtomicSpecies(
                    ID=next_id,
                    symbol=atomic_symbol,
                    aliases=atomic_aliases[(atomic_symbol, phase, state)],
                    phase=phase,
                    state=state,
                )
            )
            next_id += 1
        for template_name, phase, state in sorted(template_aliases):
            species.append(
                MolecularTemplate(
                    name=template_name,
                    aliases=template_aliases[(template_name, phase, state)],
                    phase=phase,
                    state=state,
                    matching_molecules=[],
                )
            )
        for template_name, phase in sorted(mineral_template_aliases):
            species.append(
                MineralTemplate(
                    ID=next_id,
                    name=template_name,
                    aliases=mineral_template_aliases[(template_name, phase)],
                    phase=phase,
                    template_type="wildcard_surface",
                    matching_minerals=[],
                )
            )
            next_id += 1
        for mineral_name, phase in sorted(mineral_aliases):
            species.append(
                MineralSpecies(
                    ID=next_id,
                    name=mineral_name,
                    aliases=mineral_aliases[(mineral_name, phase)],
                    phase=phase,
                )
            )
            next_id += 1
        for molecule_name, phase, state, conformation in sorted(molecule_aliases):
            species.append(
                ReferenceMolecule(
                    ID=next_id,
                    sequence=molecule_name,
                    aliases=molecule_aliases[(molecule_name, phase, state, conformation)],
                    phase=phase,
                    state=state,
                    conformation=conformation,
                )
            )
            next_id += 1
        _set_template_matching_molecules(species)
        return cls(species=species)

    @property
    def molecules(self):
        return [
            species
            for species in self.species
            if isinstance(species, (ReferenceMolecule, AtomicSpecies))
        ]
    @property
    def templates(self):
        return [
            species
            for species in self.species
            if isinstance(species, MolecularTemplate)
        ]
    @property
    def minerals(self):
        return [
            species
            for species in self.species
            if isinstance(species, MineralSpecies)
        ]

    def molecule_names(self):
        return [
            molecule.show_sequence()
            if hasattr(molecule, "show_sequence")
            else molecule.symbol
            for molecule in self.molecules
        ]

    def template_names(self):
        return [template.name for template in self.templates]

    def mineral_names(self):
        return [mineral.name for mineral in self.minerals]

#
#  define if species is molecular template
#

def _is_molecular_template(species_name):
    return (
        "_n" in species_name
        or "nplus1" in species_name
        or species_name.startswith("R_")
        or _is_lumped_polymer_template(species_name)
        or species_name in {
            "aminonitrile_precursor",
            "C_rich_polymer",
            "hydrocarbon_modified_tholin",
            "mixed_organic_residue",
            "mixed_tholin",
            "N_PAH",
            "N_rich_polymer",
            "N_rich_soluble_polymer",
            "nitrile_modified_tholin",
            "PAH_large",
            "PAH_seed",
            "prebiotic_residue",
            "refractory_hydrocarbon_residue",
            "refractory_tholin",
            "soluble_organics",
            "SOOT",
            "tholin",
        }
        or "*" in species_name
    )

#
#  define if species -> mineral
#

def _is_mineral_template(species_name):
    return species_name == "mineral*"

def _is_mineral_species(species_name):
    return is_mineral(species_name)


def _canonical_formula(species_name):
    if species_name in NAMED_MOLECULE_SEQUENCES:
        return NAMED_MOLECULE_SEQUENCES[species_name]
    tokens = re.findall(r"([A-Z][a-z]?|D)(\d*)", species_name)
    if not tokens:
        return species_name
    element_symbols = {
        element.symbol
        for element in elements
        if getattr(element, "symbol", None)
    }
    if any(symbol != "D" and symbol not in element_symbols for symbol, _ in tokens):
        return species_name
    parsed = "".join(symbol + count for symbol, count in tokens)
    if parsed != species_name:
        return species_name

    canonical_tokens = []
    for symbol, count_text in tokens:
        count = int(count_text) if count_text else 1
        if canonical_tokens and canonical_tokens[-1][0] == symbol:
            canonical_tokens[-1][1] += count
        else:
            canonical_tokens.append([symbol, count])

    canonical = ""
    for symbol, count in canonical_tokens:
        canonical += symbol
        if count > 1:
            canonical += str(count)
    return canonical


def _resolve_max_polymer_length(max_polymer_length, molecule_aliases):
    if max_polymer_length is not None:
        return int(max_polymer_length)
    max_detected_length = 0
    for molecule_name, _, _, _ in molecule_aliases:
        match = re.fullmatch(r"P\[[^\]]+\](\d+)", molecule_name)
        if match is not None:
            max_detected_length = max(max_detected_length, int(match.group(1)))
    return max_detected_length


def _add_template_polymer_molecules(template_aliases, molecule_aliases, max_polymer_length):
    if max_polymer_length < 2:
        return
    for template_name, phase, state in template_aliases:
        for monomer in TEMPLATE_POLYMER_FAMILIES.get(template_name, ()):
            for polymer_length in range(2, max_polymer_length + 1):
                molecule_key = (
                    f"P[{_canonical_formula(monomer)}]{polymer_length}",
                    phase,
                    state,
                    "linear",
                )
                molecule_aliases.setdefault(molecule_key, set())
        for sequence, conformation in PAH_TEMPLATE_REPRESENTATIVES.get(template_name, ()):
            molecule_key = (sequence, phase, state, conformation)
            molecule_aliases.setdefault(molecule_key, set())


def _set_template_matching_molecules(species):
    molecules = [
        species_item
        for species_item in species
        if isinstance(species_item, ReferenceMolecule)
    ]
    templates = [
        species_item
        for species_item in species
        if isinstance(species_item, MolecularTemplate)
    ]
    for template in templates:
        monomers = TEMPLATE_POLYMER_FAMILIES.get(template.name)
        matching_sequences = set()
        if monomers is not None:
            matching_sequences.update(
                molecule.sequence
                for molecule in molecules
                if _matches_polymer_family(molecule.sequence, monomers)
            )
        if template.name in PAH_TEMPLATE_NAMES:
            matching_sequences.update(
                (molecule.sequence, molecule.conformation)
                for molecule in molecules
                if _matches_pah_template(template.name, molecule.sequence, molecule.conformation)
            )
        if not matching_sequences:
            continue
        template.matching_molecules = [
            molecule
            for molecule in molecules
            if (
                (
                    molecule.sequence in matching_sequences
                    or (molecule.sequence, molecule.conformation) in matching_sequences
                )
                and molecule.phase == template.phase
                and molecule.state == template.state
            )
        ]


def _matches_polymer_family(sequence, monomers):
    for monomer in monomers:
        if re.fullmatch(rf"P\[{re.escape(_canonical_formula(monomer))}\]\d+", sequence):
            return True
    return False


def _matches_pah_template(template_name, sequence, conformation):
    return (sequence, conformation) in PAH_TEMPLATE_REPRESENTATIVES.get(template_name, ())


def _is_lumped_polymer_template(species_name):
    if species_name.startswith("P_CxHy_"):
        return True
    if not species_name.startswith("P_"):
        return False
    if re.search(r"_(long|frag|n|nplus1)$", species_name):
        return True
    return not re.search(r"_\d+$", species_name)


def _is_numbered_polymer_species(species_name):
    return re.fullmatch(r"P_[A-Za-z0-9]+_\d+", species_name) is not None


def _canonical_polymer_formula(species_name):
    match = re.fullmatch(r"P_([A-Za-z0-9]+)_(\d+)", species_name)
    if match is None:
        return species_name
    monomer, count = match.groups()
    return f"P[{_canonical_formula(monomer)}]{count}"


def _molecular_template_type(species_name):
    if species_name == "M":
        return "third_body"
    if species_name in {"hnu", "hv", "photon"}:
        return "radiation"
    if "_n" in species_name or "nplus1" in species_name:
        return "polymer_chain"
    if species_name.startswith("R_"):
        return "functional_group"
    if "*" in species_name:
        return "wildcard_surface"
    return "unknown"
