import argparse
# initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-ct", nargs=1, default=["evol"])
parser.add_argument("-json_input", nargs=1, default=None)
parser.add_argument("-json_NN_input", nargs=1, default=None)

# plot fitness data
parser.add_argument("-acf_data", nargs='+', default=None)
parser.add_argument("-acf_data_out", nargs='+', default=None)

# plot fitness distrib.
parser.add_argument("-size_acfs", nargs='+', default=None)
parser.add_argument("-size_cat", nargs='+', default=None)

# aij
parser.add_argument("-aij_data", nargs=1, default=None)
parser.add_argument("-aij_out", nargs=1, default=None)

# mu(t)
parser.add_argument("-mut_inp", nargs=1, default=None)
parser.add_argument("-mut_out", nargs=1, default=None)

# fitness compare
parser.add_argument("-X_size", nargs='+', default=None)

# x(t) plots
parser.add_argument("-xt_data", nargs='+', default=None)
parser.add_argument("-xt_data_out", nargs='+', default=None)