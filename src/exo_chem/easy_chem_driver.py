import numpy as np
from dataclasses import dataclass
import easychem.easychem as ec
from src.common.units import Q_
from src.utilities.logging_module import log

@dataclass
class EasyChemResult:
    mole_fractions: dict[str, float]
    pressure: Q_
    temperature: Q_
    # get abundance
    def get_mole_fraction(self, species: str, default=None):
        return self.mole_fractions.get(species, default)
    def get_abundance(self, species: str, default=None):
        return self.get_mole_fraction(species, default)
    def get_partial_pressure(self, species: str, default=None):
        mole_fraction = self.get_mole_fraction(species)
        if mole_fraction is None:
            return default
        return mole_fraction * self.pressure
    # for number density use ideal gas law

@dataclass
class EasyChemProfileResult:
    species: list[str]
    layers: list[EasyChemResult]
    mole_fraction_profiles: dict[str, np.ndarray]

#
#  check layer is valid
#

def _is_valid_easychem_result(result: EasyChemResult, chemical_species: list[str]) -> bool:
    for species in chemical_species:
        mole_fraction = result.get_mole_fraction(species, 0.0)
        if not np.isfinite(mole_fraction) or mole_fraction < 0.0:
            return False
    return True

#
#   interpolate failed layer profile
#

def _interpolate_failed_profile_layers(
        pressure,
        layers: list[EasyChemResult],
        mole_fraction_profiles: dict[str, np.ndarray],
        valid_layers: np.ndarray,
    ) -> None:
    failed_layers = ~valid_layers
    if not np.any(failed_layers):
        return
    if not np.any(valid_layers):
        log.error("EasyChem failed for every atmospheric layer; cannot interpolate chemistry profile")
    pressure_grid = np.asarray(pressure.to("bar").magnitude, dtype=float)
    interpolation_grid = np.log(pressure_grid)
    valid_grid = interpolation_grid[valid_layers]
    sort_idx = np.argsort(valid_grid)
    valid_grid = valid_grid[sort_idx]
    failed_indices = np.flatnonzero(failed_layers)
    log.warning(
        "EasyChem returned invalid mole fractions for "
        f"{failed_indices.size} layer(s); interpolating those layers in log-pressure space"
    )
    for species, profile in mole_fraction_profiles.items():
        valid_values = np.asarray(profile[valid_layers], dtype=float)[sort_idx]
        if valid_values.size == 1:
            repaired_values = np.full(failed_indices.size, valid_values[0], dtype=float)
        else:
            repaired_values = np.interp(
                interpolation_grid[failed_layers],
                valid_grid,
                valid_values,
            )
        profile[failed_layers] = np.clip(repaired_values, 0.0, 1.0)
        for ilayer in failed_indices:
            layers[ilayer].mole_fractions[species] = profile[ilayer]

#
#   find suspect invalid layers
#

def _mark_suspicious_profile_layers(
        mole_fraction_profiles: dict[str, np.ndarray],
        valid_layers: np.ndarray,
        absolute_jump_tol=5.0e-2,
        relative_jump_fraction=0.5,
    ) -> np.ndarray:
    suspicious_layers = np.zeros_like(valid_layers, dtype=bool)
    n_layers = valid_layers.size
    species_profiles = list(mole_fraction_profiles.values())
    if species_profiles:
        total_profile = np.sum(np.vstack(species_profiles), axis=0)
        profiles_to_check = species_profiles + [total_profile]
    else:
        profiles_to_check = []
    for ilayer in range(1, n_layers - 1):
        if not (valid_layers[ilayer - 1] and valid_layers[ilayer] and valid_layers[ilayer + 1]):
            continue
        for profile in profiles_to_check:
            left = profile[ilayer - 1]
            center = profile[ilayer]
            right = profile[ilayer + 1]
            neighbor_scale = max(abs(left), abs(right), 1.0e-30)
            neighbor_min = min(left, right)
            neighbor_max = max(left, right)
            outside_neighbor_range = (
                center < neighbor_min - absolute_jump_tol
                or center > neighbor_max + absolute_jump_tol
            )
            relative_jump = (
                abs(center - 0.5*(left + right))
                > max(absolute_jump_tol, relative_jump_fraction*neighbor_scale)
            )
            if outside_neighbor_range and relative_jump:
                suspicious_layers[ilayer] = True
                break
    return suspicious_layers

#
#   renormalize mole fractions -> \sum_i x_i = 1
#

def _normalize_mole_fraction_profiles(
        mole_fraction_profiles: dict[str, np.ndarray],
        layers: list[EasyChemResult],
    ) -> None:
    species_profiles = list(mole_fraction_profiles.values())
    if not species_profiles:
        return
    total_profile = np.sum(np.vstack(species_profiles), axis=0)
    bad_layers = total_profile <= 0.0
    if np.any(bad_layers):
        log.warning(
            "Total mole fraction sums to zero for "
            f"{np.flatnonzero(bad_layers).size} layer(s); "
            "normalization skipped for those layer(s)"
        )
    normalization = np.ones_like(total_profile)
    valid_layers = ~bad_layers
    normalization[valid_layers] = 1.0/total_profile[valid_layers]
    for species, profile in mole_fraction_profiles.items():
        profile *= normalization
        for ilayer, value in enumerate(profile):
            if not bad_layers[ilayer]:
                layers[ilayer].mole_fractions[species] = value

#
#   smoothen mole fractions
#

def _smooth_mole_fraction_profiles(
        mole_fraction_profiles: dict[str, np.ndarray],
        window: int = 3,
    ) -> None:
    if window < 3:
        return
    kernel = np.ones(window, dtype=float) / window
    for species, profile in mole_fraction_profiles.items():
        if profile.size < window:
            continue
        smoothed = np.convolve(profile, kernel, mode="same")
        profile[:] = np.clip(smoothed, 0.0, 1.0)

#
#   run easy chem backend
#

def run_easychem_backend(
        pressure: Q_,
        temperature: Q_,
        atomic_abund=None,
        chemical_species=None,
        log_selected=True
    ) -> EasyChemResult:
    """
    Run EasyChem once for the given thermodynamic state.
    Return species abundances.
    """
    # -------------------------------
    # Initialize solver
    # -------------------------------
    exo = ec.ExoAtmos()
    # -------------------------------
    # Set elemental abundances
    # -------------------------------
    if log_selected:
        log.info("\t Default atoms: " + str(exo.atoms))
    abunds = np.array(
        [atomic_abund.get(atom, 0.0) for atom in exo.atoms],
        dtype=float
    )
    exo.updateAtomAbunds(abunds)
    # set env. factors
    P = pressure.to("bar").magnitude
    T = temperature.to("K").magnitude
    # -------------------------------
    # Solve equilibrium
    # -------------------------------
    exo.solve(P, T)
    # -------------------------------
    # Get results
    # -------------------------------
    abundances = dict(exo.result_mol())
    if chemical_species is not None:
        abundances = {
            species: abundances.get(species, 0.0)
            for species in chemical_species
        }
    if log_selected:
        selected_species = chemical_species or ["H2", "H2O", "CO", "CH4", "NH3", "HCN", "O2"]
        log.info("\n\n")
        log.info("\t --- Selected species ---")
        for sp in selected_species:
            if sp in abundances:
                log.info(f"\t {sp:4s} : " + str(abundances[sp]))
            else:
                log.warning(f"{sp:4s} : not present")
    # return result
    return EasyChemResult(
        mole_fractions=abundances,
        pressure=pressure,
        temperature=temperature
    )

#
#   run easy chem over an altitude profile
#

def run_easy_chem_full_profile(
        pressure,
        temperature,
        atomic_abund=None,
        chemical_species=None
    ) -> EasyChemProfileResult:
    # chemical species
    if chemical_species is None:
        chemical_species = ["H2", "H2O", "CO", "CH4", "NH3", "HCN", "O2"]
    chemical_species = list(chemical_species)
    # set up output variables
    layers = []
    mole_fraction_profiles = {
        species: np.full_like(temperature.magnitude, np.nan, dtype=float)
        for species in chemical_species
    }
    valid_layers = np.zeros_like(temperature.magnitude, dtype=bool)
    # run over altitude layers
    for ilayer, (p, t) in enumerate(zip(pressure, temperature)):
        # get easy_chem result
        result = run_easychem_backend(
            pressure=p,
            temperature=t,
            atomic_abund=atomic_abund,
            chemical_species=chemical_species,
            log_selected=False,
        )
        layers.append(result)
        valid_layers[ilayer] = _is_valid_easychem_result(result, chemical_species)
        for species in chemical_species:
            mole_fraction_profiles[species][ilayer] = result.get_mole_fraction(species, 0.0)
    # get suspicious layers to recompute
    suspicious_layers = _mark_suspicious_profile_layers(
        mole_fraction_profiles=mole_fraction_profiles,
        valid_layers=valid_layers,
    )
    if np.any(suspicious_layers):
        log.warning(
            "EasyChem returned isolated suspicious profile jump(s) for "
            f"{np.flatnonzero(suspicious_layers).size} layer(s); interpolating those layers"
        )
        valid_layers[suspicious_layers] = False
    # 1) interpolate between layers
    _interpolate_failed_profile_layers(
        pressure=pressure,
        layers=layers,
        mole_fraction_profiles=mole_fraction_profiles,
        valid_layers=valid_layers,
    )
    # 2) smooth mole fractions
    _smooth_mole_fraction_profiles(
        mole_fraction_profiles=mole_fraction_profiles,
        window=3,
    )
    _normalize_mole_fraction_profiles(
        mole_fraction_profiles=mole_fraction_profiles,
        layers=layers,
    )
    # return easy chem results
    return EasyChemProfileResult(
        species=chemical_species,
        layers=layers,
        mole_fraction_profiles=mole_fraction_profiles
    )
