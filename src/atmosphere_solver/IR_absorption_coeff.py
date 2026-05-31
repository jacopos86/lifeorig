import contextlib
import io
import numpy as np
from hapi import db_begin, fetch, absorptionCoefficient_Voigt
from src.common.units import Q_
from src.utilities.plot_absorption_cross_section import plot_absorption_cross_section
from src.utilities.logging_module import log

HITRAN_SPECIES = {
    "CO2": {
        "molecule_id": 2,
        "isotope_id": 1,
        "wavenumber_ranges": [(Q_(500.0, "1 / cm"), Q_(8000.0, "1 / cm"))],
        "wavenumber_step": Q_(10.0, "1 / cm"),
    },
    "CO": {
        "molecule_id": 5,
        "isotope_id": 1,
        "wavenumber_ranges": [
            (Q_(1800.0, "1 / cm"), Q_(2300.0, "1 / cm")),
            (Q_(4000.0, "1 / cm"), Q_(4400.0, "1 / cm")),
        ],
        "wavenumber_step": Q_(5.0, "1 / cm"),
    },
    "H2O": {
        "molecule_id": 1,
        "isotope_id": 1,
        "wavenumber_ranges": [
            (Q_(500.0, "1 / cm"), Q_(12000.0, "1 / cm")),
            (Q_(12000.0, "1 / cm"), Q_(25000.0, "1 / cm")),
        ],
        "wavenumber_step": Q_(10.0, "1 / cm"),
    },
    "CH4": {
        "molecule_id": 6,
        "isotope_id": 1,
        "wavenumber_ranges": [
            (Q_(1000.0, "1 / cm"), Q_(9500.0, "1 / cm")),
            (Q_(11000.0, "1 / cm"), Q_(14000.0, "1 / cm")),
        ],
        "wavenumber_step": Q_(10.0, "1 / cm"),
    },
    "O2": {
        "molecule_id": 7,
        "isotope_id": 1,
        "wavenumber_ranges": [
            (Q_(12900.0, "1 / cm"), Q_(13200.0, "1 / cm")),
            (Q_(14300.0, "1 / cm"), Q_(14600.0, "1 / cm")),
        ],
        "wavenumber_step": Q_(1.0, "1 / cm"),
    },
    "NH3": {
        "molecule_id": 11,
        "isotope_id": 1,
        "wavenumber_ranges": [(Q_(500.0, "1 / cm"), Q_(8000.0, "1 / cm"))],
        "wavenumber_step": Q_(10.0, "1 / cm"),
    },
    "OH": {
        "molecule_id": 13,
        "isotope_id": 1,
        "wavenumber_ranges": [
            (Q_(2500.0, "1 / cm"), Q_(4000.0, "1 / cm")),
            (Q_(13000.0, "1 / cm"), Q_(16000.0, "1 / cm")),
        ],
        "wavenumber_step": Q_(5.0, "1 / cm"),
    },
    "NO": {
        "molecule_id": 8,
        "isotope_id": 1,
        "wavenumber_ranges": [
            (Q_(1700.0, "1 / cm"), Q_(2000.0, "1 / cm")),
            (Q_(5200.0, "1 / cm"), Q_(5700.0, "1 / cm")),
        ],
        "wavenumber_step": Q_(5.0, "1 / cm"),
    },
    "NO2": {
        "molecule_id": 10,
        "isotope_id": 1,
        "wavenumber_ranges": [(Q_(10000.0, "1 / cm"), Q_(25000.0, "1 / cm"))],
        "wavenumber_step": Q_(10.0, "1 / cm"),
    },
    "N2O": {
        "molecule_id": 4,
        "isotope_id": 1,
        "wavenumber_ranges": [(Q_(500.0, "1 / cm"), Q_(5000.0, "1 / cm"))],
        "wavenumber_step": Q_(10.0, "1 / cm"),
    },
    "HCN": {
        "molecule_id": 23,
        "isotope_id": 1,
        "wavenumber_ranges": [(Q_(500.0, "1 / cm"), Q_(8000.0, "1 / cm"))],
        "wavenumber_step": Q_(10.0, "1 / cm"),
    },
}

#
#    IR HITRAN/HAPI line absorption cross section sigma_i,lambda(T,P)
#
#    This module is for thermal-IR molecular line absorption. It is not the
#    shortwave/UV-visible opacity model; stellar-band attenuation should use
#    Rayleigh scattering plus UV/VIS cross-section data.
#

def hitran_species_absorption_cross_section(species, wavelength_grid, temperature, pressure):
    temperature_K = np.asarray(temperature.to("kelvin").magnitude, dtype=float)
    if species not in HITRAN_SPECIES:
        return np.zeros(
            (wavelength_grid.magnitude.size, temperature_K.size),
            dtype=float
        )
    try:
        db_begin("hitran_data")
    except Exception:
        return None
    # Convert the shared wavelength grid to HITRAN's wavenumber convention.
    wavelength_cm = wavelength_grid.to("cm").magnitude
    wavenumber = 1.0/wavelength_cm
    wavenumber_min = float(np.min(wavenumber))
    wavenumber_max = float(np.max(wavenumber))
    if wavenumber_min >= 10000.0:
        log.warning(
            "HITRAN line absorption called on a shortwave grid "
            f"({wavenumber_min:.1f}-{wavenumber_max:.1f} cm^-1); "
            "returning zero because this module is for IR opacity"
        )
        return np.zeros(
            (wavenumber.size, temperature_K.size),
            dtype=float
        )
    # Keep only the hardcoded HITRAN line ranges that overlap this spectrum.
    hitran_data = HITRAN_SPECIES[species]
    hitran_ranges = []
    for range_min, range_max in hitran_data["wavenumber_ranges"]:
        range_min = range_min.to("1 / cm").magnitude
        range_max = range_max.to("1 / cm").magnitude
        overlap_min = max(wavenumber_min, range_min)
        overlap_max = min(wavenumber_max, range_max)
        if overlap_min < overlap_max:
            hitran_ranges.append((overlap_min, overlap_max))
    if not hitran_ranges:
        return np.zeros(
            (wavenumber.size, temperature_K.size),
            dtype=float
        )
    wavenumber_step = hitran_data["wavenumber_step"].to("1 / cm").magnitude
    molecule_id = hitran_data["molecule_id"]
    isotope_id = hitran_data["isotope_id"]
    pressure_atm = np.asarray(pressure.to("atm").magnitude, dtype=float)
    sigma = np.zeros((wavenumber.size, temperature_K.size), dtype=float)
    # Run HAPI on each high-resolution molecular range, then project back to
    # the shared stellar wavelength grid.
    for irange, (range_min, range_max) in enumerate(hitran_ranges):
        table_name = f"{species}_{irange}"
        try:
            fetch(table_name, molecule_id, isotope_id, range_min, range_max)
        except Exception as exc:
            log.warning(f"HITRAN fetch failed for {species} range {range_min:.3f}-{range_max:.3f} cm^-1: {exc}")
            continue
        for ilayer, (t_layer, p_layer) in enumerate(zip(temperature_K, pressure_atm)):
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    hapi_wavenumber, hapi_sigma = absorptionCoefficient_Voigt(
                        SourceTables=table_name,
                        Environment={"p": float(p_layer), "T": float(t_layer)},
                        WavenumberRange=[range_min, range_max],
                        WavenumberStep=wavenumber_step,
                        HITRAN_units=True,
                    )
            except Exception as exc:
                log.warning(f"HITRAN cross section failed for {species} layer {ilayer}: {exc}")
                continue
            # convert to m^2
            hapi_sigma_m2 = np.asarray(hapi_sigma, dtype=float)*1.0e-4
            sigma[:, ilayer] += np.interp(
                wavenumber,
                np.asarray(hapi_wavenumber, dtype=float),
                hapi_sigma_m2,
                left=0.0,
                right=0.0,
            )
    if np.any(sigma < 0.0):
        log.error(f"HITRAN cross section for {species} must be non-negative")
    return sigma

#
#    absorption coefficient k_lambda(z) = sum_i n_i(z) sigma_i,lambda(T,P)
#

def set_IR_absorption_coefficient_profile(
        wavelength_grid,
        altitude,
        species_number_density,
        temperature,
        pressure,
        output_dir=None,
    ):
    # n. wave lengths / n. z layers
    n_wavelength = wavelength_grid.magnitude.size
    n_layers = altitude.size
    active_hitran_species = {"CO2", "O2"}
    # set absorption coeff profile
    absorption_coefficient = np.zeros((n_wavelength, n_layers), dtype=float)
    for species, number_density in species_number_density.items():
        print(species, number_density)
        if species not in active_hitran_species:
            continue
        number_density = np.asarray(number_density, dtype=float)
        sigma = hitran_species_absorption_cross_section(
            species=species,
            wavelength_grid=wavelength_grid,
            temperature=temperature,
            pressure=pressure
        )
        if sigma is None:
            log.error(f"absorption cross section species:{species} not found")
        if output_dir is not None:
            plot_absorption_cross_section(
                wavelength_grid=wavelength_grid,
                sigma=sigma,
                species=species,
                output_file=f"{output_dir}/absorption_cross_section_{species}.png",
            )
        absorption_coefficient += sigma*number_density[None, :]
    return absorption_coefficient
