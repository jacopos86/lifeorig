
#
#   Metabolites class used inside reaction networks
#

class Metabolite:
    '''
    Metabolite base class
    '''
    def __init__(self, molecule, count=0, catalyst_state=None):
        self.molecule = molecule
        self.count = count
        self.catalyst_state = catalyst_state
    @property
    def is_catalyst(self):
        return self.catalyst_state is not None
    @property
    def length(self):
        return self.molecule.length
    @property
    def sequence(self):
        return self.molecule.sequence
    def clone(self):
        if self.catalyst_state is None:
            new_cstate = None
        else:
            new_cstate = self.catalyst_state.clone()
        return Metabolite(
            molecule=self.molecule,
            count=self.count,
            catalyst_state=new_cstate,
        )