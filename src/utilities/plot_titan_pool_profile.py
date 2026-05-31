import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from src.common.units import Q_
from src.environment.pool_spatial_profile import get_pool_spatial_profile_function

#
#   plot function
#

def plot_titan_pool_profile(
        output_dir: str,
        function_name: str = "middle_hill_profile",
        length=Q_(10.0, "centimeter"),
        output_filename: str = "titan_pool_spatial_profile.png",
    ) -> str:
    profile_func = get_pool_spatial_profile_function(function_name)
    x = Q_(np.linspace(0.0, length.to("centimeter").magnitude, 400), "centimeter")
    y = profile_func(x, length=length).to("millimeter")
    # figure
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.plot(x.magnitude, y.magnitude, lw=2)
    ax.fill_between(x.magnitude, 0.0, y.magnitude, alpha=0.25)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("floor elevation (mm)")
    ax.set_title(function_name)
    fig.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)