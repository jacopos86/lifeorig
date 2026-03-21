


class Protocell:
    def __init__(self, idx):
        self.id = idx
        self.X_set = None
        self.C_net = None
        self.volume = None
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