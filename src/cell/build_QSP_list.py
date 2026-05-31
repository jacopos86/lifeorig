import numpy as np
from src.cell.spher_protocell import SphProtocell

#
#   This module sets up the quasi species system
#   It will be a list of protocells
#

def set_up_empty_QSP_list(n_protocells, protocell_data):
    return ProtocellList(
        n=n_protocells, 
        protocell_data=protocell_data
    )

#
#    Protocell list class
#

class ProtocellList:
    """
    ProtocellList manages a population of protocells.
    Each protocell should be an object with:
        - X : metabolites (dict or array)
        - Y : catalysts (list of Catalyst objects)
        - V : volume
        - step(dt) : advance internal dynamics
        - divide() : handle division
    """
    def __init__(self, n, protocell_data):
        self._list = []
        # optionally initialize with n empty protocells
        for idx in range(n):
            self.add_random_protocell(
                idx=idx, 
                protocell_data=protocell_data
            )
    def add_random_protocell(self, idx, protocell_data):
        """Add a new random protocell (user-defined)"""
        # placeholder: replace with your Protocell constructor
        proto = SphProtocell(
            idx=idx, 
            radius=protocell_data.get("radius"), 
            n_shells=protocell_data.get("n_shells")
        )
        self._list.append(proto)
    def add(self, protocell):
        """Add an existing protocell object"""
        self._list.append(protocell)
    def remove(self, protocell):
        """Remove a protocell"""
        self._list.remove(protocell)
    def __len__(self):
        return len(self._list)
    def __getitem__(self, idx):
        return self._list[idx]
    def __iter__(self):
        return iter(self._list)
    # ----------------------------
    # set initial molecules set
    # ----------------------------
    def set_initial_metabolite_distr(self, X_set, metabolites_data):
        """
        Initialize molecule counts inside each protocell.
        Parameters
        ----------
        X_set : list or dict
            Unique molecule types in the simulation
        init_size : int
            Total number of molecules per protocell
        """
        distr_type = metabolites_data.get("metabolites_distr_type")
        # normalize X_set
        if isinstance(X_set, dict):
            molecules = list(X_set.values())
        else:
            molecules = list(X_set)
        # --- compute probabilities ---
        if distr_type == "uniform":
            n_types = len(molecules)
            probs = np.ones(n_types) / n_types
        elif distr_type == "length_decay":
            alpha = metabolites_data.get("decay_const")
            lengths = np.array([len(mol) for mol in molecules])
            # exponential decay
            weights = np.exp(-alpha * lengths)
            probs = weights / weights.sum()
        else:
            log.error(f"Unknown distribution mode: {distr_type}")
        return probs
    # ----------------------------
    # Population-level dynamics
    # ----------------------------
    def step_all(self, dt):
        """
        Advance all protocells by dt:
        - internal reactions
        - volume growth
        - stochastic division
        """
        new_protocells = []
        for proto in self._list:
            proto.step(dt)
            if proto.V > proto.V_divide_threshold:
                daughters = proto.divide()
                new_protocells.extend(daughters)

        # add new daughters
        self._list.extend(new_protocells)
    def transfer_metabolites(self, phi=0.1):
        """
        Perform random metabolite transfer between protocells.
        phi : fraction of donor metabolites that can be transferred
        """
        n = len(self._list)
        if n < 2:
            return  # nothing to transfer

        for i, recipient in enumerate(self._list):
            for donor_idx, donor in enumerate(self._list):
                if donor is recipient:
                    continue
                # for each metabolite
                for key in recipient.X.keys():
                    transferred = np.random.binomial(donor.X[key], phi)
                    recipient.X[key] += transferred
                    donor.X[key] -= transferred