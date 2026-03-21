


class Reaction(ABC):
    def __init__(self):
        # catalysts
        self.catalysts = None
        # kinetic parameters
        self.k = None
    @abstractmethod
    def set_initial_catalyst_rates(self):
        ''' must be implemented in subclasses'''
        raise NotImplementedError


#
#    Cleavage reaction
#

class CleavageReaction(Reaction):
    ''' Cleavage reaction class'''
    def __init__(self, r, p1, p2):
        super().__init__()
        # reactant
        self.r = r
        # products
        self.p1 = p1
        self.p2 = p2
    def set_initial_catalyst_rates(self):
        pass


#
#    ligation reaction
#

class LigationReaction(Reaction):
    ''' Ligation reaction class'''
    def __init__(self, r1, r2, p):
        super().__init__()
        # reactant
        self.r1 = r1
        self.r2 = r2
        # products
        self.p = p
    def set_initial_catalyst_rates(self):
        pass