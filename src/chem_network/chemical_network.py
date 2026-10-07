from typing import Counter
import numpy as np
from src.reactions.reaction_rate_models import build_rate_model
from src.reactions.reaction_class import ChemicalReaction, ColliderSet
from src.utilities.logging_module import log

#
#    chemical network
#

class ReactionNetwork:
    def __init__(self, species_set):
        # molecules + templates
        self.species_set = species_set
        # species to be evolved
        self.dyn_species = None
        # species -> index
        self.species_to_index = {}
        self.reaction_list = []
    def compile(self):
        self._set_dyn_species_array()
    def add(self, parsed):
        # build rate model
        rate_model = build_rate_model(parsed)
        collider_set = self._resolve_collider_set(parsed.collider)
        # add concrete reaction
        participant_template = self._find_template_prod_react(parsed)
        if participant_template is None:
           self._add_concrete_reaction(
               reaction_id=parsed.reaction_id,
               reactant_names=parsed.reactants,
               product_names=parsed.products,
               rate_model=rate_model,
               collider_set=collider_set
           )
           return
        for molecule in participant_template.matching_molecules:
            self._add_concrete_reaction(
                reaction_id=parsed.reaction_id,
                reactant_names=parsed.reactants,
                product_names=parsed.products,
                template=participant_template,
                molecule=molecule,
                rate_model=rate_model,
                collider_set=collider_set
            )
    def _find_template_prod_react(self, parsed):
        react_names = set(parsed.reactants + parsed.products)
        for template in self.species_set.templates:
            if react_names.intersection(template.aliases):
                return template
        return None
    def _resolve_collider_set(self, collider_name):
        if collider_name == None:
            return None
        for template in self.species_set.templates:
            if collider_name in template.aliases:
                indices = np.array([
                    self.dyn_species.index(molecule)
                    for molecule in template.matching_molecules
                ])
                return ColliderSet(
                    species_indices=indices,
                    efficiencies=np.ones(len(indices))
                )
        log.error(
            f"Unknown collider template: {collider_name}"
        )
    def _add_concrete_reaction(
            self, 
            reaction_id, 
            reactant_names, 
            product_names, 
            template=None, 
            molecule=None,
            rate_model=None,
            collider_set=None
        ):
        molecule_index = None
        if template is not None:
            molecule_index = self.dyn_species.index(molecule)
        reactants = Counter(
            molecule_index
            if template and name in template.aliases
            else self.species_to_index[name]
            for name in reactant_names
            if name != "hnu"
        )
        products = Counter(
            molecule_index
            if template and name in template.aliases
            else self.species_to_index[name]
            for name in product_names
            if name != "hnu"
        )
        reaction = ChemicalReaction(
            reaction_id=reaction_id,
            react_indices=np.array(
                list(reactants.keys())
            ),
            react_coeff=np.array(
                list(reactants.values())
            ),
            prod_indices=np.array(
                list(products.keys())
            ),
            prod_coeff=np.array(
                list(products.values())
            ),
            rate_model=rate_model,
            collider_set=collider_set
        )
        self.reaction_list.append(reaction)
    @property
    def num_reactions(self):
        return len(self.reaction_list)
    @property
    def num_species(self):
        return len(self.species_set)
    def _set_dyn_species_array(self):
        species_by_id = {
            species.ID: species
            for species in self.species_set.molecules
        }
        for template in self.species_set.templates:
            for species in template.matching_molecules:
                species_by_id[species.ID] = species
        self.dyn_species = tuple(
            sorted(species_by_id.values(), key=lambda species: species.ID)
        )
        self.species_to_index = {
            name: index
            for index, species in enumerate(self.dyn_species)
            for name in species.aliases
        }
        print(self.species_to_index)