from dataclasses import dataclass

#
#   Reference molecule implementation
#

@dataclass
class AtomicSpecies:
    ID: int
    symbol: str
    aliases: set[str]
    phase: str
    state: str

@dataclass
class ReferenceMolecule:
    ''' molecule loaded from a reference reaction network '''
    ID: int
    sequence: str
    aliases: set[str]
    phase: str
    state: str
    conformation: str
    def show_sequence(self):
        return self.sequence

#
#   Factory wrapper for Molecules
#

@dataclass
class MolecularTemplate:
    name: str
    template_type: str
    aliases: set[str]
    phase: str
    state: str
    matching_molecules: list[ReferenceMolecule]