import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


#
#   plot HITRAN absorption cross section
#

def plot_absorption_cross_section(
        wavelength_grid,
        sigma,
        species: str,
        output_file: str = "absorption_cross_section.png",
        max_layers: int = 5,
    ) -> str:
    wavelength_nm = wavelength_grid.to("nanometer").magnitude
    sigma = np.asarray(sigma, dtype=float)
    n_layers = sigma.shape[1]
    layer_indices = np.unique(
        np.linspace(0, n_layers - 1, min(max_layers, n_layers), dtype=int)
    )
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for ilayer in layer_indices:
        sigma_layer = sigma[:, ilayer]
        positive = sigma_layer > 0.0
        if np.any(positive):
            ax.semilogy(
                wavelength_nm[positive],
                sigma_layer[positive],
                lw=1.5,
                label=f"layer {ilayer}",
            )
        else:
            ax.plot(
                wavelength_nm,
                sigma_layer,
                lw=1.5,
                label=f"layer {ilayer}",
            )
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel(r"$\sigma_\lambda$ (m$^2$ molecule$^{-1}$)")
    ax.set_title(f"{species} absorption cross section")
    ax.grid(alpha=0.25, which="both")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(output_file, dpi=160)
    plt.close(fig)
    return output_file
