import numpy as np
import copy
from src.input_data.read_input import p

class RAFEngine:
    """
    this class
    computes maxRAF and all irrRAFs using the reaction structure.
    """
    def __init__(self, CNET):
        # set of reactions
        self.ligand_reactions = CNET.ligand_reactions
        self.cleavage_reactions = CNET.cleavage_reactions
        # n. possible network config.
        self.n_config = p.nconfig
        self.size_F = CNET.size_F
    ##################################################################
    # 1. Convert your reactions to a unified internal representation
    ##################################################################
    def _build_reaction_set(self):
        reaction_set = []
        # ligand reactions: r1 + r2 --c--> p
        for r in self.ligand_reactions:
            reaction_set.append({
                "reactants": {r["r1_int"], r["r2_int"]},
                "catalysts": r["cList"],
                "products": {r["p_int"]},
                "raw": r,
            })
        # cleavage reactions: r --c--> p1 + p2
        for r in self.cleavage_reactions:
            reaction_set.append({
                "reactants": {r["r_int"]},
                "catalysts": r["cList"],
                "products": {r["p1_int"], r["p2_int"]},
                "raw": r,
            })
        return reaction_set
    ##################################################################
    # 2. Compute closure
    ##################################################################
    def compute_closure_set(self, ic, R, F, reaction_set):
        """
        R: list of reaction indices
        F: set or list of food molecule IDs
        reaction_set: list of reaction dicts (from _build_reaction_set)
        """
        X = set(F)
        added = True
        while added:
            added = False
            for i in R:
                r = reaction_set[i]
                catal_ic = r["catalysts"][ic]
                if r["reactants"].issubset(X) and catal_ic in X:
                    new = r["products"] - X
                    if new:
                        X |= new
                        added = True
        return X
    ##################################################################
    # 3. maxRAF implementation (adapted to your structure)
    ##################################################################
    def find_maxRAF(self):
        reaction_set = self._build_reaction_set()
        for i in range(len(reaction_set)):
            print(i, reaction_set[i])
        N = len(reaction_set)
        F = set(range(1, self.size_F + 1))
        maxRafList = []
        # run over possible config.
        for ic in range(self.n_config):
            R = list(range(N))
            while True:
                old = set(R)
                Cl = self.compute_closure_set(ic, R, F, reaction_set)
                R = [
                    i for i in R
                    if reaction_set[i]["reactants"].issubset(Cl)
                    and reaction_set[i]["catalysts"][ic] in Cl
                ]
                if set(R) == old:
                    maxRafList.append(R)
                    break
        return maxRafList   # list of reaction indices
    ##################################################################
    # 4. Restricted maxRAF (needed for irrRAF search)
    ##################################################################
    def maxRAF_sub(self, R0, reaction_set, F):
        R = list(R0)

        while True:
            old = set(R)
            Cl = self.compute_closure_set(R, F, reaction_set)

            R = [
                i for i in R
                if reaction_set[i]["reactants"].issubset(Cl)
                and reaction_set[i]["catalysts"].issubset(Cl)
            ]

            if set(R) == old:
                break

        return R if R else None

    ##################################################################
    # 5. Enumerate all irreducible RAFs
    ##################################################################
    def all_irrRAFs(self):
        reaction_set = self._build_reaction_set()
        F = set(range(1, self.size_F + 1))

        max_raf = self.maxRAF()
        if not max_raf:
            return []

        return self._irrRAF_recursive(max_raf, reaction_set, F)

    def _irrRAF_recursive(self, R, reaction_set, F):
        irrRAFs = []

        for r in R:
            R_without = [x for x in R if x != r]
            sub = self.maxRAF_sub(R_without, reaction_set, F)

            if sub is not None:
                # Recursive branching: r is optional
                irrRAFs += self._irrRAF_recursive(sub, reaction_set, F)
            else:
                # r is essential — cannot be removed
                pass

        # If no smaller RAF found → R is irreducible
        if not irrRAFs:
            return [R]

        # Remove duplicates
        uniq = []
        out = []
        for raf in irrRAFs:
            key = tuple(sorted(raf))
            if key not in uniq:
                uniq.append(key)
                out.append(raf)
        return out