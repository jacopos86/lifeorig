

#
#  This module builds the catalysts set
#

def build_catalysts_set(size_X, size_C):
    catalyst_set = []
    if catal_prob_distr is not None:
        # use gaussian prob. distr.
        mean = 10  # Mean of the distribution
        std_dev = 2  # Standard deviation of the distribution
        num_samples = 100  # Number of samples to generate
        continuous_samples = np.random.normal(loc=mean, scale=std_dev, size=num_samples)
    else:
        while len(catalyst_set) < size_C:
            i = random.choice(np.arange(1, size_X+1, 1))
            if i not in catalyst_set:
                catalyst_set.append(i)
    return catalyst_set