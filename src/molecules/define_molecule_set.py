from itertools import product
from src.input_data.read_input import p
from src.molecules.molecule_model import Molecule, MultiPolymer
from src.utilities.logging_module import log

#
#   build the molecules set
#

def build_molecule_set(metabolites_data):
    """
    Build the global set of unique Molecule objects present in the simulation.
    Returns:
        molecules : list[Molecule]
        molecule_map : dict[sequence -> Molecule]
    """
    max_polymer_length = metabolites_data.get("pol_strng_maxsize")
    mol_type = metabolites_data.get("type")
    molecules = []
    molecule_map = {}
    if mol_type == "binary":
        # build all binary polymers up to max length
        for L in range(1, max_polymer_length + 1):
            for i in range(2 ** L):
                seq = format(i, f'0{L}b')  # binary string
                mol = Molecule(seq)
                molecules.append(mol)
                molecule_map[seq] = mol
    elif mol_type == 'multi':
        alphabet = list(MultiPolymer._MONOMER_TYPES)
        for L in range(1, max_polymer_length + 1):
            for seq in product(alphabet, repeat=L):
                mol = Molecule(list(seq))
                molecules.append(mol)
                molecule_map[tuple(seq)] = mol
    else:
        log.error("Unknown molecule type")
    log.info("\n")
    for i, mol in enumerate(molecules):
        log.info(f"\t molecule {i+1}: {mol.show_sequence()} -- length: {len(mol)}")
    log.info("\t " + p.sep)
    log.info("\n")
    return molecules, molecule_map