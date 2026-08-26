import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from src.common.units import Q_
from src.environment.surf_pond.pool_spatial_profile import get_pool_spatial_profile_function

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


def titan_surface_area_profile(
        z,
        height,
        surface_area,
        profile_data: dict | None = None,
    ):
    profile_data = profile_data or {"type": "uniform"}
    profile_type = profile_data.get("type", "uniform")
    z_cm = z.to("centimeter").magnitude
    height_cm = height.to("centimeter").magnitude
    surface_area_cm2 = surface_area.to("centimeter ** 2").magnitude
    if profile_type == "uniform":
        profile = np.ones_like(z_cm)
    elif profile_type == "bottom_weighted":
        decay_length = profile_data.get("decay_length")
        if decay_length is None:
            decay_length = Q_(0.1 * height_cm, "centimeter")
        ell_cm = decay_length.to("centimeter").magnitude
        profile = np.exp(-z_cm / ell_cm)
    else:
        raise ValueError(f"Unknown Titan surface area profile: {profile_type}")
    return Q_(surface_area_cm2 * profile, "centimeter ** 2")


def plot_titan_surface_area_profile(
        output_dir: str,
        height,
        surface_area,
        profile_data: dict | None = None,
        reactive_area_to_surface_area_ratio: float = 1.0,
        output_filename: str = "titan_surface_area_profile.png",
    ) -> str:
    z = Q_(np.linspace(0.0, height.to("centimeter").magnitude, 400), "centimeter")
    surface_area_profile = titan_surface_area_profile(
        z=z,
        height=height,
        surface_area=surface_area,
        profile_data=profile_data,
    )
    reactive_area_profile = surface_area_profile * reactive_area_to_surface_area_ratio
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.plot(z.magnitude, surface_area_profile.magnitude, lw=2, label="surface area")
    ax.plot(z.magnitude, reactive_area_profile.magnitude, lw=2, ls="--", label="reactive area")
    ax.set_xlabel("z (cm)")
    ax.set_ylabel("A(z) (cm^2)")
    ax.set_title((profile_data or {"type": "uniform"}).get("type", "uniform"))
    ax.legend()
    fig.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return output_path
