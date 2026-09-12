from abc import abstractmethod, ABC
from src.utilities.logging_module import log

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

#
#   Abstract base molecule class
#

class StringMoleculeModel(ABC):
    def __init__(self, sequence):
        self.sequence = sequence
        self.properties = {}  # dictionary for additional features (stability, polarity, etc.)
    @abstractmethod
    def ligate(self, other):
        """Combine with another molecule"""
        raise NotImplementedError
    @abstractmethod
    def cleave(self, index):
        """Split molecule at given index"""
        raise NotImplementedError
    @abstractmethod
    def mutate_conformation(self, mutation_rate=0.01):
        """Mutate molecule sequence"""
        raise NotImplementedError
    def compute_features(self):
        """Compute and store any features of the molecule"""
        # default: length and basic counts
        self.properties['length'] = self.length
        self.properties['composition'] = {c: self.sequence.count(c) for c in set(self.sequence)}
        self.properties['conformation'] = self.conformation
        return self.properties
    @abstractmethod
    def show_sequence(self):
        ''' return sequence as string'''
        raise NotImplementedError
    def __len__(self):
        return len(self.sequence)
    @property
    def length(self):
        return len(self.sequence)

#
#   Binary polymer implementation
#

class BinaryPolymer(StringMoleculeModel):
    ''' binary polymer class '''
    _MONOMER_TYPES = {'0', '1'}
    def __init__(self, sequence):
        super().__init__(sequence)
        # ensure only 0 and 1
        assert all(c in '01' for c in sequence), "BinaryPolymer sequence must contain only 0 or 1"
    def ligate(self, other):
        if not isinstance(other, BinaryPolymer):
            raise TypeError("Can only ligate with another BinaryPolymer")
        new_seq = self.sequence + other.sequence
        return BinaryPolymer(new_seq)
    def cleave(self, index):
        if index < 1 or index >= self.length:
            raise ValueError("Cleavage index out of range")
        return BinaryPolymer(self.sequence[:index]), BinaryPolymer(self.sequence[index:])
    def mutate_conformation(self, mutation_rate=0.01):
        pass
    def show_sequence(self):
        return self.sequence

#
#   Multi-monomer polymer (more general)
#

class MultiPolymer(StringMoleculeModel):
    ''' Multipolymer class'''
    _MONOMER_TYPES = {'A', 'B', 'C', 'D'}
    def __init__(self, sequence):
        if isinstance(sequence, str):
            # convert string to list of monomers
            sequence = str(sequence)
        super().__init__(sequence)
    def ligate(self, other):
        if not isinstance(other, MultiPolymer):
            log.error("Can only ligate with another MultiPolymer")
        return MultiPolymer(self.sequence + other.sequence)
    def cleave(self, index):
        if index < 1 or index >= self.length:
            log.error("Cleavage index out of range")
        return MultiPolymer(self.sequence[:index]), MultiPolymer(self.sequence[index:])
    def mutate_conformation(self, mutation_rate=0.01, monomer_set=None):
        pass
    def show_sequence(self):
        seq = ""
        for c in self.sequence:
            seq += str(c)
        return seq
