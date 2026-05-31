from abc import ABC, abstractmethod
from src.reactions.reaction_class import CleavageReaction, LigationReaction

#
#    abstract reaction network class
#

class ReactionNetwork(ABC):
    """
    Abstract reaction network class.

    A reaction network knows:
    - which reaction channels exist
    - how to compute their rates / propensities
    - how to apply a selected reaction to a protocell state
    """
    def __init__(self, network_type: str):
        # type of network
        self.network_type = network_type
        # max. num. reactions
        self.num_reactions = None
        # reactions set
        self.reaction_list = None
        # molecules set
        self.X_set = None
        self.species_set = None
        # X set size
        self.size_X = None
        # maps
        self.x_to_id = None
        self.id_to_x = None
        self.seq_to_id = None
        # max string len
        self.max_strng_size = None
    # --------------------------------------------------
    # interface
    # --------------------------------------------------
    def get_species_set(self):
        return self.species_set
    def get_available_reactions(self):
        return self.reaction_list
    def add_reaction(self, reaction: CleavageReaction | LigationReaction):
        """
        Add one reaction channel to the network.
        """
        self.reaction_list.append(reaction)
    def __len__(self):
        """
        Return number of reaction channels.
        """
        return self.num_reactions
    def _clear_reactions(self):
        if self.reaction_list is not None:
            self.reaction_list = []
    def get_available_reactions(self):
        """
        Return the set/list of reactions allowed.

        Returns
        -------
        list
            List of available reaction descriptors or indices.
        """
        return self.reaction_list
    @abstractmethod
    def build_network(self):
        """
        Buiild the network of reactions.
        """
        pass
    @abstractmethod
    def _validate_network(self):
        """
        Network validation procedure
        """
        pass