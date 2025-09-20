import matplotlib.pyplot as plt
import numpy as np

#
#   distributions plot functions
#

def plot_reaction_rate_distr(rr_ACFS, rr_extset, output_file):
    plt.ylim([-0.1, 1.1])
    plt.xlabel(r"catalysts distribution")
    plt.ylabel(r"reaction rates")
    xgr = rr_ACFS["xgr"]
    rr = rr_ACFS["rr"]
    plt.scatter(xgr, rr, alpha=0.5, color='red')
    xgr = rr_extset["xgr"]
    rr = rr_extset["rr"]
    plt.scatter(xgr, rr, alpha=0.5, color='red')
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.close()

#
#  plot ACFS distribution
#

def plot_ACFS_distr(frozen_rv, output_file):
    # 2. Define the range of x-values for plotting
    # It's good practice to create a range that covers most of the distribution's mass.
    # For a normal distribution, the mean ± 4 standard deviations is a common choice.
    x_values = np.linspace(frozen_rv.ppf(0.001), frozen_rv.ppf(0.999), 100)
    # 3. Calculate the PDF for each x-value using the frozen object
    pdf_values = frozen_rv.pdf(x_values)
    # 4. Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, pdf_values, label='Frozen Normal PDF', color='blue')
    plt.title('PDF of a Frozen Continuous Distribution')
    plt.xlabel('X-values')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.close()

#
#  plot ACFS distribution with histogram
#

def plot_ACFS_hist_distr(samples_int, size_X, output_file):
    plt.figure(figsize=(10, 6))
    plt.hist(samples_int, bins=np.arange(0, size_X+2)-0.5, edgecolor='black', alpha=0.7)
    plt.title('Histogram of Catalysts Distribution')
    plt.xlabel('Sample Values')
    plt.ylabel('Frequency')
    plt.xticks(np.arange(0, size_X+1))  # Set x-ticks to be integers between 0 and size_X
    plt.grid(True)
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.close()