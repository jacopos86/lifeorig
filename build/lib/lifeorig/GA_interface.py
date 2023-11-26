from lifeorig.read_input import p
from lifeorig.logging_module import log
from numpy.random import randint
from lifeorig.objective_func import objective
#
#  interface classical / quantum
#  GA
class GA_interface:
    def __init__(self):
        # n. chromosomes in population
        self.n_bits = p.n_bits
        self.n_pop = p.n_pop
    def generate_instance(self):
        from lifeorig.QGA import QGA
        from lifeorig.CGA import CGA
        if p.GA_type == "quantum":
            return QGA ()
        elif p.GA_type == "classical":
            return CGA ()
        else:
            log.error("GA algorithm: quantum or classical...")
    # set up population
    def set_initial_population(self):
        self.pop = [randint(0, 2, self.n_bits).tolist() for _ in range(self.n_pop)]
    # evaluate the score of all candidates
    # in the population
    def compute_scores(self):
        self.scores = [objective(c) for c in self.pop]