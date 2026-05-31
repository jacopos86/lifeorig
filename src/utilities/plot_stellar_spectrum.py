import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

#
#   plot function utility
#

def plot_stellar_B_lambda(
        wavelength_grid,
        B_lambda,
        output_file: str = "stellar_B_lambda.png",
    ) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.plot(
        wavelength_grid.to("nanometer").magnitude,
        B_lambda.to("W / m^3 / steradian").magnitude,
        lw=2,
    )
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel(r"$B_\lambda(T)$ (W m$^{-3}$ sr$^{-1}$)")
    ax.set_title("stellar black-body intensity")
    fig.tight_layout()
    fig.savefig(output_file, dpi=160)
    plt.close(fig)
