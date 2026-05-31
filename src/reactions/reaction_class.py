from abc import ABC, abstractmethod


class Reaction(ABC):
    def __init__(self, reaction_id=None):
        self.reaction_id = reaction_id
        # list of molecules that may catalyze the reaction
        self.catalysts = None
        # kinetic parameters
        self.k = None
    @property
    @abstractmethod
    def reactants(self):
        raise NotImplementedError
    @property
    @abstractmethod
    def products(self):
        raise NotImplementedError
    @property
    @abstractmethod
    def reaction_type(self):
        raise NotImplementedError
    def add_catalyst(self, mol, rate):
        self.catalysts.append(mol)
        self.k[mol] = rate
    def set_initial_catalyst_rates(self, catalyst_set, rate_builder):
        """
        catalyst_set  : iterable of catalyst-capable molecules
        rate_builder  : callable(rate_builder(mol, self) -> float)
        """
        raise NotImplementedError


#
#    Cleavage reaction
#

class CleavageReaction(Reaction):
    """ Cleavage reaction class """
    def __init__(self, r, p1, p2, reaction_id=None):
        super().__init__(reaction_id=reaction_id)
        # reactant
        self.r = r
        # products
        self.p1 = p1
        self.p2 = p2
    @property
    def reactants(self):
        return [self.r]
    @property
    def products(self):
        return [self.p1, self.p2]
    @property
    def reaction_type(self):
        return "cleavage"
    def set_initial_catalyst_rates(self, catalyst_set, rate_builder):
        pass


#
#    ligation reaction
#

class LigationReaction(Reaction):
    """ Ligation reaction class """
    def __init__(self, r1, r2, p, reaction_id=None):
        super().__init__(reaction_id=reaction_id)
        # reactant
        self.r1 = r1
        self.r2 = r2
        # products
        self.p = p
    @property
    def reactants(self):
        return [self.r1, self.r2]
    @property
    def products(self):
        return [self.p]
    @property
    def reaction_type(self):
        return "ligation"
    def set_initial_catalyst_rates(self, catalyst_set, rate_builder):
        pass