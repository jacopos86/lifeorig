from numpy.random import randint
#
#  tournament selection
def tournament_selection(chromosomes, k=3):
        # first random selection
        selection_ix = randint(chromosomes.n_pop)
        for ix in randint(0, chromosomes.n_pop, k-1):
            if chromosomes.scores[ix] < chromosomes.scores[selection_ix]:
                selection_ix = ix
        return chromosomes.pop[selection_ix]