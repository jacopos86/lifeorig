from abc import ABC, abstractmethod
import numpy as np
from src.common.units import Q_, internal_units
from src.common.phys_constants import G
from src.utilities.logging_module import log
from src.utilities.plot_chemical_mass_profiles import plot_chemical_mass_profiles
from src.exo_chem.easy_chem_driver import run_easy_chem_full_profile
from src.atmosphere_solver.hydrostatic_solver import set_hydro_solver
from src.atmosphere_solver.radiative_solver import RadiativeSolver
from src.atmosphere_solver.IR_absorption_coeff import set_IR_absorption_coefficient_profile
from src.atmosphere_solver.optical_depth_eval import (
    OpticalDepthProfiles,
    absorption_optical_depth_profile,
)
from src.atmosphere_solver.atmosph_data import AtmLayerDyn

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
        T0 = self.atmosphere_data.get("initial_temperature")
        self._initial_temperature = (
            288.0
            if T0 is None
            else float(T0.to(self._units.get("temperature")).magnitude)
        )
        # set hydrostatic solver
        self._hydro_solver = set_hydro_solver(
            atmosphere_data,
            planet_data,
            units=self._units,
            atomic_abundances=self.chem_input.atomic_abundances,
            chemical_species=self.chem_input.chemical_species
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
    # build initial state
    def _build_initial_state_layers(self):
        # this is very Earth like
        # altitude
        z = self._build_altitude_grid()
        # temperature
        temperature = self._initial_temperature_profile(z)
        # density
        density = self._initial_density_profile(z)
        # gravity profile
        gravity = self._gravity_profile(z)
        # set surface pressure
        self._set_surface_pressure()
        # pressure evaluation
        pressure = self._solve_hydrostatic_pressure(
            z=z,
            temperature=temperature,
            density=density,
            gravity=gravity
        )
        # compute chemical composition
        chem_data = run_easy_chem_full_profile(
            pressure=Q_(pressure, self.pressure_unit),
            temperature=Q_(temperature, self.temperature_unit),
            atomic_abund=self.chem_input.atomic_abundances,
            chemical_species=self.chem_input.chemical_species
        )
        # molecular weights
        mu = self._mean_molecular_weight_profile(chem_data, z)
        
        # species number density
        species_number_density = self._species_number_density(
            density=density,
            chem_data=chem_data,
            mu=mu
        )
        return AtmLayerDyn(
            altitude=z,
            pressure=pressure,
            temperature=temperature,
            mean_molecular_mass=mu,
            gravity=gravity,
            species_number_density=species_number_density,
        )
    # build altitude grid
    def _build_altitude_grid(self):
        n_layers = int(self.atmosphere_data["n_layers"])
        z_min = 0.0
        z_max = self.atmosphere_data["z_max"].to(self._units["length"]).magnitude
        # Internal solver variables are unitless floats in the units set by
        # _set_internal_units().
        return np.linspace(
            z_min,
            z_max,
            n_layers,
        )
    # gravity profile
    def _gravity_profile(self, z):
        planet_mass = self.planet_data.planet_mass.to(self._units["mass"]).magnitude
        planet_radius = self.planet_data.planet_radius.to(self._units["length"]).magnitude
        grav_const = G.to(f"{self._units['length']}^3 / {self._units['mass']} / {self._units['time']}^2").magnitude
        return grav_const*planet_mass/(planet_radius + z)**2
    # set initial temperature profile
    def _initial_temperature_profile(self, z):
        temperature_surface_guess_K = self._initial_temperature
        return np.full_like(z, temperature_surface_guess_K, dtype=float)
    #
    #    validation section
    #
    def _validate_input(self):
        if self.chem_input.atomic_abundances is None:
            log.error("atomic_abundances is required for layered atmosphere solvers")
        if self.stellar_data is None:
            log.error("stellar_params is required for layered atmosphere solvers")
        if self.planet_data is None:
            log.error("planetary_data is required for layered atmosphere solvers")
        if self.atmosphere_data is None:
            log.error("planetary_data.atmosphere is required for layered atmosphere solvers")
        if "n_layers" not in self.atmosphere_data:
            log.error("planetary_data.atmosphere.n_layers is required")
        if "z_max" not in self.atmosphere_data:
            log.error("planetary_data.atmosphere.z_max is required")
        # n. layers
        n_layers = self.atmosphere_data["n_layers"]
        if not isinstance(n_layers, int) or n_layers < 2:
            log.error("planetary_data.atmosphere.n_layers must be an integer >= 2")
        # z max altitude
        z_max = self.atmosphere_data["z_max"]
        if z_max is None or not hasattr(z_max, "to"):
            log.error("planetary_data.atmosphere.z_max must be a quantity")
        if z_max <= Q_(0.0, "m"):
            log.error("planetary_data.atmosphere.z_max must be > 0")
        settings = self._hydro_solver.settings
        if settings.max_iter < 1:
            log.error("hydrostatic_solver.settings.max_iter must be >= 1")
        if settings.abs_tol < 0.0:
            log.error("hydrostatic_solver.settings.abs_tol must be non-negative")
        if settings.rel_tol < 0.0:
            log.error("hydrostatic_solver.settings.rel_tol must be non-negative")
        if settings.logp_tol < 0.0:
            log.error("hydrostatic_solver.settings.logp_tol must be non-negative")
        if settings.damping <= 0.0 or settings.damping > 1.0:
            log.error("hydrostatic_solver.settings.damping must satisfy 0 < damping <= 1")
        if settings.min_pressure <= 0.0:
            log.error("hydrostatic_solver.settings.min_pressure must be > 0")
        if settings.anderson_depth < 1:
            log.error("hydrostatic_solver.settings.anderson_depth must be >= 1")
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
