

# =========================================================
#   single catalyst class
# =========================================================

class Catalyst:
    def __init__(self, n_reactions, 
                 init_mean=0.0, 
                 init_std=1.0,
                 mutation_std=0.1):
        """
        n_reactions : number of reactions catalyzed
        init_mean   : mean of log-rate initial distribution
        init_std    : std of log-rate initial distribution
        mutation_std: std of Gaussian mutation in log space
        """
        self.n_reactions = n_reactions
        self.mutation_std = mutation_std
        # work in log space
        self.log_rates = np.random.normal(
            loc=init_mean,
            scale=init_std,
            size=n_reactions
        )
    @property
    def rates(self):
        return np.exp(self.log_rates)
    def mutate(self, mutation_rate=0.1):
        """
        Each reaction rate mutates independently
        with probability mutation_rate.
        """
        for r in range(self.n_reactions):
            if np.random.rand() < mutation_rate:
                self.log_rates[r] += np.random.normal(
                    0.0, self.mutation_std
                )
    def clone(self):
        """ Create a deep copy (for protocell division) """
        C_new = Catalyst(self.n_reactions)
        C_new.log_rates = np.copy(self.log_rates)
        C_new.mutation_std = self.mutation_std
        return C_new
    def summary(self):
        log.info("\t " + p.sep)
        log.info("\n")
        log.info("\t catalyst rate distribution : ")
        log.info("\n")
        line = ""
        for j in range(self.n_reactions):
            line += " {0:.3f}".format(np.exp(self.log_rates[j]))
        log.info("\t " + line)
        log.info("\n")
        log.info("\t " + p.sep)