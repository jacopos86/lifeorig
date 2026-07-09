from itertools import product
from src.input_data.chemical_network_parser import ChemicalNetworkParser
from src.metabolites.metabolite_class import Metabolite
from src.chemical_types.molecule_model import Molecule
from src.chemical_types.molecule_model import MultiPolymer
from src.input_data.chemical_species_parser import parse_species_into_chemicals
from src.utilities.logging_module import log

#
#  1) build metabolites driver function
#

def build_metabolites(metabolites_params):
    mol_type = metabolites_params.get("type")
    if mol_type == "binary":
        return build_binary_metabolites(metabolites_params)
    if mol_type == "multi":
        return build_multi_metabolites(metabolites_params)
    if mol_type == "reference_file":
        return build_reference_file_metabolites(metabolites_params)
    log.error(f"Unknown molecule type: {mol_type}")

#
#  2) set binary metabolites model
#

def build_binary_metabolites(metabolites_params):
    max_size = metabolites_params.get("pol_strng_maxsize")
    molecules = []
    molecule_map = {}
    for length in range(1, max_size + 1):
        for idx in range(2 ** length):
            sequence = format(idx, f"0{length}b")
            molecule = Molecule(sequence)
            molecules.append(molecule)
            molecule_map[sequence] = molecule
    initial_population = build_uniform_initial_population(
        molecules,
        metabolites_params.get("initial_population_molecules", 0),
    )
    return molecules, molecule_map, initial_population

#
#  3) SET multi mode metabolites model
#

def build_multi_metabolites(metabolites_params):
    max_size = metabolites_params.get("pol_strng_maxsize")
    molecules = []
    molecule_map = {}
    alphabet = list(MultiPolymer._MONOMER_TYPES)
    for length in range(1, max_size + 1):
        for sequence in product(alphabet, repeat=length):
            molecule = Molecule(list(sequence))
            molecules.append(molecule)
            molecule_map[tuple(sequence)] = molecule
    initial_population = build_uniform_initial_population(
        molecules,
        metabolites_params.get("initial_population_molecules", 0),
    )
    return molecules, molecule_map, initial_population

#
#  4)  build metabolites list from reference file
#

def build_reference_file_metabolites(metabolites_params):
    reaction_file = metabolites_params.get("reaction_file")
    network_data = ChemicalNetworkParser(reaction_file).parse()
    log.info(f"\t species list: {network_data.species}")
    species_names = network_data.species
    chemical_set = parse_species_into_chemicals(species_names)
    print(chemical_set)
    exit()
    molecule_map = {molecule.show_sequence(): molecule for molecule in molecules}
    initial_population = build_uniform_initial_population(
        molecules,
        metabolites_params.get("initial_population_molecules", 0),
    )
    return molecules, molecule_map, initial_population

#
#   5)   build the initial population of metabolites
#

def build_uniform_initial_population(molecules, total_count):
    if not molecules:
        return {}
    base_count = total_count // len(molecules)
    remainder = total_count % len(molecules)
    population = {}
    for idx, molecule in enumerate(molecules):
        count = base_count + (1 if idx < remainder else 0)
        population[molecule] = Metabolite(molecule=molecule, count=count)
    return population
