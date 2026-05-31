import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from src.common.units import Q_


#
#   plot atmospheric pressure profile
#

def plot_pressure_profile(
        pressure,
        altitude=None,
        pressure_unit: str = "bar",
        length_unit: str = "km",
        output_file: str = "pressure_profile.png",
    ) -> str:
    p = Q_(np.asarray(pressure, dtype=float), "Pa").to(pressure_unit).magnitude
    n_layers = p.size
    if altitude is None:
        x = np.arange(n_layers)
        x_label = "layer"
    else:
        x = Q_(np.asarray(altitude, dtype=float), "m").to(length_unit).magnitude
        x_label = f"altitude ({length_unit})"
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.semilogy(x, p, color="tab:blue", lw=2.0)
    ax.set_xlabel(x_label)
    ax.set_ylabel(f"pressure ({pressure_unit})")
    ax.set_title("pressure profile")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(output_file, dpi=160)
    plt.close(fig)
    return output_file
