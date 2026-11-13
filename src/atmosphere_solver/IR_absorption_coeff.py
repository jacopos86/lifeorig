import contextlib
import hashlib
import io
import json
import numpy as np
import os
import tempfile
from pathlib import Path
from hapi import db_begin, fetch, absorptionCoefficient_Voigt
from src.parallelization.mpi import mpi
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
    wavenumber = 1.0 / wavelength_cm
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
#   build IR cross section table
#

def build_IR_cross_section_table(
        *,
        species,
        temperature_grid,
        pressure_grid,
        wavelength_grid,
        broadening,
        spectroscopy
):
    """Return sigma[T_index, P_index, wavelength_index] in m²/molecule.

    broadening: HAPI diluent fractions, e.g. {"air": 1.0}.
    spectroscopy: numerical settings, including wavenumber_step_cm1.
    Called only on MPI root by get_IR_cross_section_table().
    """
    temperatures = np.asarray(
        temperature_grid.to("K").magnitude, dtype=float
    )
    pressures = np.asarray(
        pressure_grid.to("atm").magnitude, dtype=float
    )
    wavelengths_cm = np.asarray(
        wavelength_grid.to("cm").magnitude, dtype=float
    )
    # HAPI evaluates spectra in wavenumber space.
    target_wavenumbers = 1.0 / wavelengths_cm
    # check species data
    species_data = HITRAN_SPECIES.get(species)
    if species_data is None:
        return np.zeros(
            (
                temperatures.size,
                pressures.size,
                target_wavenumbers.size
            ),
            dtype=float
        )
    # wavelengths
    nu_min = float(target_wavenumbers.min())
    nu_max = float(target_wavenumbers.max())
    step = float(spectroscopy["wavenumber_step_cm1"])
    wing = float(spectroscopy.get("wing_cm1", 25.0))
    # Restrict this initial implementation to the standard broadening data.
    if not broadening or set(broadening) - {"air", "self"}:
        log.error(
            "This implementation supports air/self broadening only"
        )
    if any(value < 0.0 for value in broadening.values()):
        log.error("Broadening fractions must be nonnegative")
    if not np.isclose(sum(broadening.values()), 1.0):
        log.error("Broadening fractions must sum to one")
    print("OK HERE", species)
    sigma = np.zeros(
        (temperatures.size, pressures.size, target_wavenumbers.size),
        dtype=float
    )
    # Select overlapping molecular spectral ranges.
    hitran_ranges = []
    for range_min, range_max in species_data["wavenumber_ranges"]:
        overlap_min = max(nu_min, range_min.to("1 / cm").magnitude)
        overlap_max = min(nu_max, range_max.to("1 / cm").magnitude)
        if overlap_min < overlap_max:
            hitran_ranges.append((overlap_min, overlap_max))
    if not hitran_ranges:
        return sigma
    print("ok", hitran_ranges)
    with contextlib.redirect_stdout(io.StringIO()):
        db_begin("hitran_data")
    molecule_id = species_data["molecule_id"]
    isotope_id = species_data["isotope_id"]
    # run over wave numbers
    for irange, (range_min, range_max) in enumerate(hitran_ranges):
        table_name = f"{species}_{irange}"
        # Download transitions once for this spectral interval.
        with contextlib.redirect_stdout(io.StringIO()):
            fetch(
                table_name,
                molecule_id,
                isotope_id,
                max(0.0, range_min - wing),
                range_max + wing
            )
        # Reuse those transitions for every temperature–pressure pair.
        for iT, temperature in enumerate(temperatures):
            log.info(
                f"Building {species}: T={temperature:g} K "
                f"({iT + 1}/{temperatures.size})"
            )
            for iP, pressure in enumerate(pressures):
                with contextlib.redirect_stdout(io.StringIO()):
                    hapi_wavenumber, hapi_sigma = absorptionCoefficient_Voigt(
                        SourceTables=table_name,
                        Environment={
                            "T": float(temperature),
                            "p": float(pressure)
                        },
                        Diluent=broadening,
                        WavenumberRange=[range_min, range_max],
                        WavenumberStep=step,
                        WavenumberWing=wing,
                        WavenumberWingHW=0.0,
                        HITRAN_units=True
                    )
                # Convert cm²/molecule to m²/molecule.
                hapi_sigma_m2 = (
                    np.asarray(hapi_sigma, dtype=float) * 1.0e-4
                )
                # Store on the requested wavelength grid.
                sigma[iT, iP, :] += np.interp(
                    target_wavenumbers,
                    np.asarray(hapi_wavenumber, dtype=float),
                    hapi_sigma_m2,
                    left=0.0,
                    right=0.0
                )
    finite = np.isfinite(sigma)
    print("Non-finite values:", np.count_nonzero(~finite))
    print("Negative values:", np.count_nonzero(sigma < 0.0))
    #print("Finite min/max:", sigma[finite].min(), sigma[finite].max())
    bad = np.argwhere(~finite | (sigma < 0.0))
    if bad.size:
        iT, iP, iw = bad[0]
        print(
            "First invalid value:",
            "T =", temperatures[iT],
            "p_atm =", pressures[iP],
            "wavenumber =", target_wavenumbers[iw],
            "sigma =", sigma[iT, iP, iw],
        )
    #if not np.all(np.isfinite(sigma)) or np.any(sigma < 0.0):
    #    log.error(f"Invalid IR cross sections for {species}")
    return sigma

#
#    get IR cross section table
#

def get_IR_cross_section_table(
        *,
        species,
        temperature_grid,
        pressure_grid,
        wavelength_grid,
        broadening,
        spectroscopy,
        cache_dir
):
    """Collectively load or build one species' IR cross-section table.

    Grids are Pint quantities.
    broadening and spectroscopy are JSON-serializable configuration dicts.
    Returns NumPy cross sections in m²/molecule, shaped (n_T, n_P, n_lambda).
    All MPI ranks must call this with the same arguments.
    """
    temperatures = np.asarray(
        temperature_grid.to("K").magnitude, dtype=float
    )
    pressures = np.asarray(
        pressure_grid.to("Pa").magnitude, dtype=float
    )
    wavelengths = np.asarray(
        wavelength_grid.to("nm").magnitude, dtype=float
    )
    # metadata
    metadata = {
        "schema_version": 1,
        "species": species,
        "temperature_K": temperatures.tolist(),
        "pressure_Pa": pressures.tolist(),
        "wavelength_nm": wavelengths.tolist(),
        "cross_section_units": "m^2 / molecule",
        "broadening": broadening,
        "spectroscopy": spectroscopy
    }
    metadata_json = json.dumps(
        metadata, sort_keys=True, allow_nan=False
    )
    cache_key = hashlib.sha256(
        metadata_json.encode("utf-8")
    ).hexdigest()
    # path
    path = Path(cache_dir) / f"{cache_key}.npz"
    expect_shape = (
        temperatures.size,
        pressures.size,
        wavelengths.size
    )
    # set up cross section table
    error = None
    if mpi.rank == mpi.root:
        temporary_path = None
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.exists():
                with np.load(path, allow_pickle=False) as stored:
                    if stored["metadata"].item() != metadata_json:
                        raise ValueError(f"Opacity metadata mismatch: {path}")
                    if stored["sigma"].shape != expect_shape:
                        raise ValueError(f"Opacity table shape mismatch: {path}")
                log.info(f"\t Loading IR cross-section table for {species}")
            else:
                log.info(f"\t Building IR cross-section table for {species}")
                sigma = build_IR_cross_section_table(
                    species=species,
                    temperature_grid=Q_(temperatures, "K"),
                    pressure_grid=Q_(pressures, "Pa"),
                    wavelength_grid=Q_(wavelengths, "nm"),
                    broadening=broadening,
                    spectroscopy=spectroscopy
                )
                sigma = np.asarray(sigma, dtype=float)
                if sigma.shape != expect_shape:
                    raise ValueError("Generated opacity table has wrong shape")
                #if not np.all(np.isfinite(sigma)) or np.any(sigma < 0):
                #    raise ValueError("Generated cross sections are invalid")
                # Publish only a completely written table.
                with tempfile.NamedTemporaryFile(
                    dir=path.parent, suffix=".npz", delete=False
                ) as stream:
                    temporary_path = Path(stream.name)
                    np.savez_compressed(
                        stream,
                        metadata=metadata_json,
                        sigma=sigma,
                    )
                os.replace(temporary_path, path)
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
        finally:
            if temporary_path is not None and temporary_path.exists():
                temporary_path.unlink()
    # Communicate build failures before other ranks attempt loading.
    error = mpi.comm.bcast(error, root=mpi.root)
    if error is not None:
        log.error(f"IR table preparation failed: {error}")
    # Requires a cache directory accessible to every rank.
    sigma = None
    error = None
    try:
        with np.load(path, allow_pickle=False) as stored:
            sigma = stored["sigma"].copy()
    except Exception as exc:
        error = f"Rank {mpi.rank}: {type(exc).__name__}: {exc}"
    errors = mpi.comm.allgather(error)
    if any(message is not None for message in errors):
        raise RuntimeError(f"IR table loading failed: {errors}")
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
