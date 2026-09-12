import os
from src.utilities.logging_module import log
from src.input_data.chemical_network_parser import ChemicalNetworkParser
from src.input_data.vpl_atmos_reaction_parser import VPLAtmosReactionParser
from src.utilities.file_hash import _file_sha256
from src.network_generation.reaction_mysql_db import delete_old_reaction_file_data, open_reaction_database, reaction_file_is_current, insert_reaction_file_data
from src.global_parameters.global_variables import PROJECT_DIR

REFERENCE_REACTIONS_DIR = PROJECT_DIR / "reference_reactions"

#
#   main reaction DB driver
#

def reaction_database_driver(input_params):
    """
    Prepare the reaction database for the current run.

    This does not solve chemistry.
    It only makes sure the selected reaction networks are available
    in the reaction database before the molecular species set is built.
    """
    reaction_source_files = _select_reaction_files(input_params)
    # update reactions
    for reaction_file in reaction_source_files:
        _register_update_reaction_file(reaction_file)
    return reaction_source_files

#
#   select the reaction files
#

def _select_reaction_files(input_params):
    if input_params.planetary_data.name.lower() == "titan":
        titan_reaction_dir = REFERENCE_REACTIONS_DIR / "TITAN"
        return [
            titan_reaction_dir / "titan_gas_phase_hebrard2013.txt",
            titan_reaction_dir / "titan_photolysis_selected.txt",
            titan_reaction_dir / "titan_liquid_surface_polymerization.txt",
        ]
    reaction_data = input_params.chemical_network_data
    if reaction_data and reaction_data.get("reaction_files"):
        return [
            REFERENCE_REACTIONS_DIR / reaction_file
            for reaction_file in reaction_data["reaction_files"]
        ]
    if reaction_data and reaction_data.get("reaction_file"):
        return [REFERENCE_REACTIONS_DIR / reaction_data["reaction_file"]]
    return []

#
#   register reaction data
#

def _register_update_reaction_file(reaction_file):
    log.info("\t REGISTERING REACTION DATA: " + str(reaction_file))
    if not reaction_file.exists():
        raise FileNotFoundError(f"Reaction source file does not exist: {reaction_file}")
    file_hash = _file_sha256(reaction_file)
    # database
    db = open_reaction_database()
    try:
        force_import = os.environ.get("LIFEORIG_FORCE_REACTION_DB_IMPORT") == "1"
        if not force_import and reaction_file_is_current(db, reaction_file, file_hash):
            log.info("\t reaction data already current")
            return
        if reaction_file.name == "reactions.rx":
            parsed_network = VPLAtmosReactionParser(reaction_file).parse()
        else:
            parsed_network = ChemicalNetworkParser(reaction_file).parse()
        delete_old_reaction_file_data(db, reaction_file)
        insert_reaction_file_data(db, reaction_file, file_hash, parsed_network)
        db.commit()
        log.info("\t imported reactions: " + str(len(parsed_network.reactions)))
        log.info("\t imported species: " + str(len(parsed_network.species)))
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
