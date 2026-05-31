


#  class describing the
#  reaction network -> we use a multi polymer model

class MultiReactionNetwork(ReactionNetwork):
    def __init__(self):
        """
        Multi polymer model reaction network class
        """
        super().__init__(
            network_type="binary_polymer"
        )