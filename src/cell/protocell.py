import numpy as np
from abc import ABC, abstractmethod
from src.molecules.metabolite_class import Metabolite
from src.catalysts.catalyst_class import CatalystState
from src.chem_network.bin_reaction_network import BinaryReactionNetwork

#
#   Protocell class
#

class Protocell(ABC):
    """
    Protocell class
    """
    def __init__(self, idx, contact_to_mineral_surface):
        self.id = idx
        # metabolites data
        self.M_set = None
        self.X_set_init_popul = None
        # reactions network
        self.Rnet = None
        # protocell structural data
        self.Rcm = None
        self.volume = None
        self.Np = None
        self.mineral_interact = contact_to_mineral_surface
    def set_catalytic_state(self, X_set, Y_set):
        """
        Create dictionary of metabolites for this protocell
        with catalytic state

        X_set : list of all molecule types
        Y_set : subset of X_set that can act as catalysts
        """
        self.M_set = {}
        for mol in X_set:
            # assign catalyst state only if in Y_set
            if mol in Y_set:
                catalyst_state = CatalystState(
                    n_reactions=n_reactions,
                    init_mean=init_mean,
                    init_std=init_std,
                    mutation_std=mutation_std,
                )
            else:
                catalyst_state = None
            # define metabolite
            self.M_set[mol] = Metabolite(
                molecule=mol,
                count=self.X_set_init_popul[mol],
                catalyst_state=catalyst_state,
            )
        log.info("\t metabolites initialization complete")
        log.info("\t " + p.sep)
    def set_initial_metabolite_population(self, metabolites_data, X_set, metabol_distr):
        """
        Set up initial metabolite population
        """
        self.Np = metabolites_data.get("initial_population_molecules")
        # random multinomial distribution
        counts = np.random.multinomial(self.Np, metabol_distr)
        assert len(counts) == len(X_set)
        assert sum(counts) == self.Np
        # set up initial population
        self.X_set_init_popul = {
            mol: count
            for mol, count in zip(X_set, counts)
                if count > 0
        }
    def set_initial_protocell_state(self, rng=None):
        """
        Create the list of particle objects from the initial metabolite population.
        Each particle has a metabolite type and a relative position.
        This defines together with protocell radius and center of mass and volume
        the initial protocell state
        """
        if rng is None:
            rng = np.random.default_rng()
        # check whether initial X_set population is set
        if self.X_set_init_popul is None:
            log.error("Initial metabolite population must be set first.")
        # initialize particle list
        self.particles = []
        # sample relative positions in sphere
        for mol, count in self.X_set_init_popul.items():
            coord_list = self.sample_relative_metabolite_positions(count, rng=rng)
            for xyz in coord_list:
                particle = Particle(
                    metabolite_type=mol,
                    Rp=xyz
                )
                self.particles.append(particle)
        assert len(self.particles) == self.Np
    @abstractmethod
    def sample_relative_metabolite_positions(count, rng):
        raise NotImplementedError("sample_relative_metabolite_positions")
    def set_reaction_network(self, metabolites_params, X_set):
        """
        Set individual reaction network
        """
        if metabolites_params.get("type") == "binary":
            self.Rnet = BinaryReactionNetwork()
        elif metabolites_params.get("type") == "multi":
            self.Rnet = MultiReactionNetwork()
        else:
            log.error("Wrong type metabolic network")
        # build reaction network
        self.Rnet.build_network(X_set)
        self.Rnet.summary()
    def divide(self, mutation_matrix):
        daughter1 = deepcopy(self)
        daughter2 = deepcopy(self)
        # Binomial partition metabolites
        for i in range(len(self.molecules)):
            x = self.molecules[i]
            x1 = np.random.binomial(x, 0.5)
            daughter1.molecules[i] = x1
            daughter2.molecules[i] = x - x1
        # Optional mutation of catalytic maps
        #mutate_catalytic_rates(daughter1.catalytic_network, rate=mu)
        #mutate_catalysis(daughter2.catalytic_network, rate=mu)
        return daughter1, daughter2
    def info(self):
        """
        print protocell info stats
        """
        log.info("\n\t " + p.sep)
        log.info(f"\t PROTOCELL ID: {self.id}")
        log.info(f"\t PROTOCELL VOLUME: {self.volume}")
        log.info(f"\t PROTOCELL COM: {self.Rc}")
        log.info(f"\t PROTOCELL RADIUS: {self.radius}")
        log.info(f"\t NUMBER PARTICLES: {self.Np}")