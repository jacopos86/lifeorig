import os
from src.atmosphere_solver.atm_struct_driver import AtmosphSolver
from src.network_generation.reaction_mysql_db import (
    get_reaction_network_from_database,
    open_reaction_database,
)
from src.network_generation.reaction_database_driver import (
    REFERENCE_REACTIONS_DIR,
    _register_update_reaction_file,
)
from src.input_data.chemical_network_parser import ParsedChemicalNetwork
from src.utilities.plot_chemical_mass_profiles import plot_chemical_mass_profiles
from src.chemical_types.define_molecule_set import MolecularSpeciesSet
from src.utilities.logging_module import log
from src.chem_network.chemical_network import ReactionNetwork
from src.reactions.reaction_class import ChemicalReaction

#
#   Non-equilibrium layered structure
#

class NonEquilLayeredAtmosphSolver(AtmosphSolver):
    def __init__(self, chem_input, stellar_data, planet_data, atmosphere_data, output_dir=None):
        super().__init__(
            chem_input=chem_input,
            stellar_data=stellar_data,
            planet_data=planet_data,
            atmosphere_data=atmosphere_data,
            output_dir=output_dir,
        )
        self._chemical_network = None
        self._Rnet = None
        # set atmospheric chemical network
        self._set_chemical_network()
        self._set_reaction_network()
        exit()
    # run main driver
    def run(self):
        # 1) set base arrays
        z = self._build_altitude_grid()
        g = self._gravity_profile(z)
        # 2) set initial arrays guess
        T0_z = self._initial_temperature_profile(z)
        # 3) set initial variables
        layered_variables = self._hydro_solver.initialize(
            altitude=z,
            temperature=T0_z,
            gravity=g
        )
        plot_chemical_mass_profiles(
            mu=layered_variables.mean_molecular_mass,
            chem_data=layered_variables.chemistry,
            altitude=layered_variables.altitude,
            mass_unit=self._units.get("mass"),
            length_unit="km",
            output_file=os.path.join(self.output_dir, "chemical_mass_profiles_initial_guess.png"),
        )
    # set chemical network
    def _set_chemical_network(self):
        reaction_files = self.chem_input.reaction_source_files
        if not reaction_files:
            raise ValueError(
                "layered_disequilibrium requires exo_chemistry.reaction_files"
            )
        reaction_source_files = [
            REFERENCE_REACTIONS_DIR / reaction_file
            for reaction_file in reaction_files
        ]
        for reaction_file in reaction_source_files:
            _register_update_reaction_file(reaction_file)
        db = open_reaction_database()
        try:
            self._chemical_network = self._merge_reaction_networks(
                db, reaction_source_files
            )
        finally:
            db.close()
    @staticmethod
    def _merge_reaction_networks(db, reaction_source_files):
        """Merge sources in input order, retaining only the first duplicate."""
        reactions = []
        species = []
        seen_reactions = set()
        seen_species = set()
        # run over reaction files
        for reaction_file in reaction_source_files:
            network = get_reaction_network_from_database(
                db, source_files=[reaction_file]
            )
            for reaction in network.reactions:
                reaction_key = (
                    tuple(sorted(reaction.reactants)),
                    tuple(sorted(reaction.products)),
                    reaction.reversible,
                )
                if reaction_key in seen_reactions:
                    continue
                seen_reactions.add(reaction_key)
                reactions.append(reaction)
                for species_name in reaction.reactants + reaction.products:
                    if species_name in {"hnu", "M"}:
                        continue
                    if species_name not in seen_species:
                        species.append(species_name)
                        seen_species.add(species_name)
        return ParsedChemicalNetwork(species=species, reactions=reactions)
    # reaction network
    def _set_reaction_network(self):
        chemical_species = self._set_chemical_species()
        self._Rnet = ReactionNetwork(
            species=chemical_species.molecules
        )
        for parsed_reaction in self._chemical_network.reactions:
            reaction = ChemicalReaction(
                reaction_id=parsed_reaction.reaction_id,
                reactants=parsed_reaction.reactants,
                products=parsed_reaction.products,
                reaction_type=parsed_reaction.catalyst_or_control
            )
            self._Rnet.add(reaction)
            exit()
