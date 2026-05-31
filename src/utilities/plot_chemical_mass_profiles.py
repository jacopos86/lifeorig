import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from periodictable import formula
from src.common.units import Q_
from src.utilities.logging_module import log


def _smooth_profile(profile: np.ndarray, window: int = 5) -> np.ndarray:
    if window < 3 or profile.size < window:
        return profile
    pad = window // 2
    padded = np.pad(profile, pad_width=pad, mode="reflect")
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(padded, kernel, mode="valid")

#
#   plot chemical mass contribution profiles
#

def plot_chemical_mass_profiles(
        mu,
        chem_data,
        altitude=None,
        mass_unit: str = "kg",
        length_unit: str = "km",
        output_file: str = "chemical_mass_profiles.png",
        max_species: int = 8,
        smooth_window: int = 5,
    ) -> str:
    mu_amu = Q_(np.asarray(mu, dtype=float), mass_unit).to("amu").magnitude
    mu_amu = _smooth_profile(mu_amu, window=smooth_window)
    n_layers = mu_amu.size
    if altitude is None:
        x = np.arange(n_layers)
        x_label = "layer"
    else:
        x = Q_(np.asarray(altitude, dtype=float), "m").to(length_unit).magnitude
        x_label = f"altitude ({length_unit})"
    contributions = {}
    for species, mole_fraction in chem_data.mole_fraction_profiles.items():
        try:
            molecular_weight_amu = formula(species).mass
        except Exception as exc:
            log.warning(f"cannot plot mass contribution for species {species}: {exc}")
            continue
        xi = np.asarray(mole_fraction, dtype=float)
        contribution = xi*molecular_weight_amu
        contribution = _smooth_profile(contribution, window=smooth_window)
        if np.any(np.isfinite(contribution)):
            contributions[species] = contribution
    species_to_plot = sorted(
        contributions,
        key=lambda sp: np.nanmax(contributions[sp]),
        reverse=True,
    )[:max_species]
    fig, (ax_mu, ax_species) = plt.subplots(
        2,
        1,
        figsize=(8.0, 6.0),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.6]},
    )
    ax_mu.plot(x, mu_amu, color="black", lw=2.4)
    ax_mu.set_ylabel(r"$\mu$ (amu)")
    ax_mu.set_title("chemical mass profile")
    ax_mu.grid(alpha=0.25)
    for species in species_to_plot:
        ax_species.plot(x, contributions[species], lw=1.6, label=species)
    ax_species.set_xlabel(x_label)
    ax_species.set_ylabel(r"$x_i m_i$ (amu)")
    ax_species.grid(alpha=0.25)
    ax_species.legend(loc="best", fontsize=8, ncol=2)
    fig.tight_layout()
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(output_file, dpi=160)
    plt.close(fig)
    return output_file

def plot_chemical_mass_profile_lines(
        mu,
        chem_data,
        altitude=None,
        mass_unit: str = "kg",
        length_unit: str = "km",
        output_file: str = "chemical_mass_profile_lines.png",
        max_species: int = 8,
        smooth_window: int = 5,
    ) -> str:
    mu_amu = Q_(np.asarray(mu, dtype=float), mass_unit).to("amu").magnitude
    mu_amu = _smooth_profile(mu_amu, window=smooth_window)
    n_layers = mu_amu.size
    if altitude is None:
        x = np.arange(n_layers)
        x_label = "layer"
    else:
        x = Q_(np.asarray(altitude, dtype=float), "m").to(length_unit).magnitude
        x_label = f"altitude ({length_unit})"
    contributions = {}
    for species, mole_fraction in chem_data.mole_fraction_profiles.items():
        try:
            molecular_weight_amu = formula(species).mass
        except Exception as exc:
            log.warning(f"cannot plot mass contribution for species {species}: {exc}")
            continue
        xi = np.asarray(mole_fraction, dtype=float)
        contribution = xi*molecular_weight_amu
        contribution = _smooth_profile(contribution, window=smooth_window)
        if np.any(np.isfinite(contribution)):
            contributions[species] = contribution
    species_to_plot = sorted(
        contributions,
        key=lambda sp: np.nanmax(contributions[sp]),
        reverse=True,
    )[:max_species]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot(x, mu_amu, color="black", lw=2.4, label=r"$\mu$")
    for species in species_to_plot:
        ax.plot(x, contributions[species], lw=1.5, label=species)
    ax.set_xlabel(x_label)
    ax.set_ylabel("mass contribution (amu)")
    ax.set_title("chemical mass profile")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8, ncol=2)
    fig.tight_layout()
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(output_file, dpi=160)
    plt.close(fig)
    return output_file
