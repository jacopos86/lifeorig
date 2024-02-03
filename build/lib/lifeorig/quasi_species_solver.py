# The quasi species solver solves
# the quasi species equation given in input
# 1) a list of fitness values
# 2) a random matrix Q describing the mutations
class QuasiSpeciesSolver():
    def __init__(self, Q, f):
        # random matrix 
        self.Q = Q
        # fitness
        self.f = f