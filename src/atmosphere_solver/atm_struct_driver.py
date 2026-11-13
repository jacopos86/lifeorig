from abc import ABC, abstractmethod
import numpy as np
from pathlib import Path
from src.parallelization.array_backend import to_host
from src.common.units import Q_, internal_units, get_magnitude
from src.utilities.logging_module import log
from src.utilities.plot_chemical_mass_profiles import plot_chemical_mass_profiles
from src.exo_chem.easy_chem_driver import run_easy_chem_full_profile
from src.atmosphere_solver.hydrostatic_solver import set_hydro_solver, hydrostatic_driver
from src.atmosphere_solver.radiative_solver import RadiativeSolver
from src.atmosphere_solver.IR_absorption_coeff import set_IR_absorption_coefficient_profile
from src.atmosphere_solver.optical_depth_eval import (
    OpticalDepthProfiles,
    absorption_optical_depth_profile,
)
from src.atmosphere_solver.initial_atmospheric_state import InitialAtmosphericStateBuilder
from src.atmosphere_solver.atmosph_data import AtmLayerDyn, IROpacitySettings

#
#   full atmospheric solver
#

class AtmosphSolver(ABC):
    def __init__(self, chem_input, stellar_data, planet_data, atmosphere_data, output_dir=None):
        self.chem_input = chem_input
        self.planet_data = planet_data
        self.stellar_data = stellar_data
        self.atmosphere_data = atmosphere_data
        self.output_dir = output_dir or "."
        # set internal units
        self._set_internal_units()
        # set atmosph. mass
        self._atmosph_mass = None
        self._set_atmosph_mass()
        # Initial grid, temperature, gravity and chemical composition.
        self._inital_state_builder = InitialAtmosphericStateBuilder(
            chem_input=self.chem_input,
            planet_data=self.planet_data,
            atmosph_data=self.atmosphere_data,
            units=self._units
        )
        # set hydrostatic solver
        self._hydro_solver = set_hydro_solver(
            atmosphere_data,
            units=self._units
        )
        # here set radiative solver
        self._radiative_solver = RadiativeSolver(
            atmosphere_data,
            stellar_data,
            planet_data,
            units=self._units,
            output_dir=self.output_dir
        )
        # validate data
        self._validate_input()
    # set internal units
    def _set_internal_units(self):
        self._units = internal_units
    # atmosph. mass
    def _set_atmosph_mass(self):
        atmosph_mass_frac = float(self.atmosphere_data.get("atmosphere_mass_fraction"))
        self._atmosph_mass = get_magnitude(
            value=atmosph_mass_frac * self.planet_data.planet_mass,
            units=self._units.get("mass"),
            default_value=None
        )
    # build initial state
    def _build_initial_state_layers(self, *, pressure_guess=None):
        """Construct initial composition and consistent hydrostatic layers."""
        method = self.chem_input.initial_chemistry["method"]
        equilibrium_update = None
        # easychem method
        if method == "easychem":
            if pressure_guess is None:
                log.error(
                    "EasyChem initialization requires a pressure profile guess "
                    "in internal pressure units"
                )
            # Called only by the initial chemistry–hydrostatic consistency loop.
            def equilibrium_update(*, pressure, temperature):
                chemistry = run_easy_chem_full_profile(
                    pressure=Q_(
                        to_host(pressure),
                        self._units["pressure"]
                    ),
                    temperature=Q_(
                        to_host(temperature),
                        self._units["temperature"]
                    ),
                    atomic_abund=self.chem_input.atomic_abundances,
                    chemical_species=self.chem_input.chemical_species
                )
                return chemistry.mole_fraction_profiles
        # 1) get initial configuration
        # Profile: load composition and temperature.
        # Prescribed: construct composition on the configured grid.
        # EasyChem: calculate the first composition using the pressure guess.
        initial_state = self._inital_state_builder.build(
            network=self.chem_system.reaction_network,
            pressure_guess=pressure_guess
        )
        # prepare opacity tables
        ir_opacity = IROpacitySettings.from_input(
            self.atmosphere_data["opacity_tables"]["IR"]
        )
        self._radiative_solver.prepare_opacity_tables(
            initial_state["species"],
            temperature_grid=ir_opacity.temperature_grid,
            pressure_grid=ir_opacity.pressure_grid,
            broadening=ir_opacity.broadening,
            spectroscopy=ir_opacity.spectroscopy,
            cache_dir=(Path(self.output_dir) / ir_opacity.cache_dir).resolve(),
            absorbers=ir_opacity.absorbers
        )
        # 2) solve hydrostatic eqs.
        # return completed initial layers
        self._hydro_solver.set_chem_species(initial_state["species"])
        return hydrostatic_driver(
            solver=self._hydro_solver,
            state=initial_state,
            atmosph_mass=self._atmosph_mass,
            planet_radius=float(
                self.planet_data.planet_radius.to(
                    self._units["length"]
                ).magnitude
            ),
            initialize=True,
            initialization_method=method,
            equilibrium_update=equilibrium_update
        )
    #
    #    validation section
    #
    def _validate_input(self):
        """Validate common solver inputs and hydrostatic settings."""
        if self.chem_input is None:
            raise ValueError("chem_input is required")
        if self.stellar_data is None:
            raise ValueError("stellar_params is required for layered atmosphere solvers")
        if self.planet_data is None:
            raise ValueError("planetary_data is required for layered atmosphere solvers")
        if not isinstance(self.atmosphere_data, dict):
            raise ValueError("planetary_data.atmosphere must be a mapping")
        settings = self._hydro_solver.settings
        if settings.max_iter < 1:
            raise ValueError("hydrostatic_solver.settings.max_iter must be >= 1")
        for name in ("abs_tol", "rel_tol", "logp_tol"):
            value = getattr(settings, name)
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"hydrostatic_solver.settings.{name} must be finite and non-negative")
        if not np.isfinite(settings.damping) or not 0.0 < settings.damping <= 1.0:
            raise ValueError("hydrostatic_solver.settings.damping must satisfy 0 < damping <= 1")
        if not np.isfinite(settings.min_pressure) or settings.min_pressure <= 0.0:
            raise ValueError("hydrostatic_solver.settings.min_pressure must be finite and > 0")
        if settings.anderson_depth < 1:
            raise ValueError("hydrostatic_solver.settings.anderson_depth must be >= 1")
    # check T / P profiles
    def _check_temperature_pressure_profiles(self, temperature, pressure):
        if np.any(temperature <= 0.0):
            log.error("temperature must be positive for equilibrium atmosphere solve")
        if np.any(pressure <= 0.0):
            log.error("pressure must be positive for equilibrium atmosphere solve")
    #
    #   MAIN DRIVER RUN
    #
    @abstractmethod
    def run(self):
        raise NotImplementedError
