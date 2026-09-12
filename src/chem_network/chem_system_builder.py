from dataclasses import dataclass
from src.chemical_types.molecule_set_builder import build_molecular_species_set
from src.chem_network.chemical_network import ReactionNetwork
from src.network_generation.reaction_database_driver import (
    REFERENCE_REACTIONS_DIR,
    _register_update_reaction_file,
)
from src.network_generation.reaction_mysql_db import (
    get_reaction_network_from_database,
    open_reaction_database,
)
from src.utilities.logging_module import log

@dataclass
class ChemicalSystem:
    species_set: object
    reaction_network: ReactionNetwork

#
#   build chemical system
#

def build_chemical_system(chem_input, input_params):
    reaction_files = [
        REFERENCE_REACTIONS_DIR / filename
        for filename in chem_input.reaction_source_files
    ]

    for reaction_file in reaction_files:
        _register_update_reaction_file(reaction_file)

    db = open_reaction_database()
    try:
        parsed_network = get_reaction_network_from_database(
            db,
            source_files=reaction_files,
        )
    finally:
        db.close()

    #
    #   build species set
    #

    species_set = build_molecular_species_set(
        input_params,
        reaction_source_files=reaction_files,
    )

    log.info("Molecules:")
    for molecule in species_set.molecules:
        log.info("\n")
        log.info(f"\t{molecule}")

    log.info("Templates:")
    for template in species_set.templates:
        log.info("\n")
        log.info(f"\t{template}")

    log.info("Minerals:")
    for mineral in species_set.minerals:
        log.info("\n")
        log.info(f"\t{mineral}")

    #
    # build reaction network
    #
     
    reaction_network = ReactionNetwork(species_set)
    reaction_network.compile()

    for parsed_reaction in parsed_network.reactions:
        print(parsed_reaction)
        reaction_network.add(parsed_reaction)
        print(reaction_network.reaction_list)
        exit()