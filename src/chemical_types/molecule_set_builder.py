from src.chemical_types.define_molecule_set import MolecularSpeciesSet
from src.network_generation.reaction_mysql_db import (
    get_species_from_reaction_database,
    open_reaction_database,
)

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
