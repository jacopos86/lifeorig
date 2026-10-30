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
from src.parallelization.mpi import mpi

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
    registration_error = None
    if mpi.rank == mpi.root:
        try:
            for reaction_file in reaction_files:
                _register_update_reaction_file(reaction_file)
        except Exception as error:
            registration_error = (
                type(error).__name__,
                str(error),
            )
    registration_error = mpi.comm.bcast(
        registration_error,
        root=mpi.root,
    )
    if registration_error is not None:
        error_type, error_message = registration_error
        log.error(
            "Reaction database registration failed on "
            f"MPI rank {mpi.root}: "
            f"{error_type}: {error_message}"
        )
    mpi.comm.Barrier()
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
    # RUN OVER REACTIONS HERE
    unsupported_reactions = []
    for parsed_reaction in parsed_network.reactions:
        try:
            start = reaction_network.num_reactions
            reaction_network.add(parsed_reaction)
            generated = reaction_network.num_reactions - start
            if generated > 1:
                participant_names = set(
                    parsed_reaction.reactants + parsed_reaction.products
                )
                matched_templates = [
                    template.name
                    for template in species_set.templates
                    if participant_names.intersection(template.aliases)
                ]
                log.info(
                    f"Expanded {parsed_reaction.reaction_id}: "
                    f"{parsed_reaction.equation} using templates "
                    f"{matched_templates} -> {generated} channels"
                )
        except NotImplementedError as error:
            unsupported_reactions.append(
                (parsed_reaction.reaction_id, str(error))
            )
    if unsupported_reactions:
        log.warning(
            f"{len(unsupported_reactions)} reactions are not implemented"
        )
        for reaction_id, error_msg in unsupported_reactions:
            log.warning(f"{reaction_id}: {error_msg}")
    log.info(
        f"Built {reaction_network.num_reactions} reaction channels "
        f"from {len(parsed_network.reactions)} parsed reactions"
    )
    return ChemicalSystem(
        species_set=species_set,
        reaction_network=reaction_network
    )