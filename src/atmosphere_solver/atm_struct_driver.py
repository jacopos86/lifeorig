from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
import numpy as np
from periodictable import formula
from scipy.integrate import cumulative_trapezoid
from src.common.units import Q_
from src.common.phys_constants import G, kb
from src.utilities.logging_module import log
from src.utilities.plot_chemical_mass_profiles import plot_chemical_mass_profiles
from src.utilities.plot_stellar_spectrum import plot_stellar_B_lambda
from src.exo_chem.easy_chem_driver import run_easy_chem_full_profile
from src.atmosphere_solver.IR_absorption_coeff import set_IR_absorption_coefficient_profile
from src.atmosphere_solver.optical_depth_eval import (
    OpticalDepthProfiles,
    absorption_optical_depth_profile,
)

@dataclass
class AtmLayerDyn:
    altitude: np.ndarray
    pressure: np.ndarray
    temperature: np.ndarray
    mean_molecular_mass: np.ndarray
    gravity: np.ndarray
    species_number_density: dict[str, np.ndarray]
    chemistry: object | None = None


@dataclass
class AtmDynResult:
    atomic_abundances: dict[str, float]
    stellar_params: object
    uv_wavelength_grid: Q_
    ir_wavelength_grid: Q_
    stellar_B_lambda: Q_
    bond_albedo: float | None
    spectral_albedo: np.ndarray | None
    layers: AtmLayerDyn

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
        # physical parameters
        self.UV_wavelength_grid = None
        self.IR_wavelength_grid = None
        self.stellar_B_lambda = None
        self.F0 = None
        self.bond_albedo = None
        self.spectral_albedo = None
        # set internal units
        self._set_internal_units()
        # loop controls
        self._max_iter_loop = int(self.atmosphere_data.get("max_iter_loop", 50))
        self._rel_tol = float(self.atmosphere_data.get("rel_tol"))
        self._abs_tol = float(self.atmosphere_data.get("abs_tol").to(self.pressure_unit).magnitude)
        self._damping_loop = float(self.atmosphere_data.get("damping_loop", 1.0))
        # set atmospheric mass
        self._atmosphere_mass_fraction = float(self.atmosphere_data.get("atmosphere_mass_fraction"))
        self._atmosphere_mass = self._atmosphere_mass_fraction*self.planet_data.planet_mass
        # pressure at surface
        self._P_surf = None
        self._min_pressure = Q_(1e-12, self.pressure_unit)
        # validate data
        self._validate_input()
    # set internal units
    def _set_internal_units(self):
        self.units = {
            "length": "m",
            "mass": "kg",
            "time": "s",
            "energy": "joule",
            "temperature": "K",
            "pressure": "Pa",
            "density": "kg / m^3",
            "number_density": "1 / m^3",
            "gravity": "m / s^2",
            "wavelength": "nanometer",
            "spectral_radiance": "W / m^3 / steradian",
            "spectral_flux": "W / m^3",
        }
        self.length_unit = self.units["length"]
        self.mass_unit = self.units["mass"]
        self.time_unit = self.units["time"]
        self.temperature_unit = self.units["temperature"]
        self.pressure_unit = self.units["pressure"]
        self.density_unit = self.units["density"]
        self.number_density_unit = self.units["number_density"]
        self.gravity_unit = self.units["gravity"]
        self._kB = kb.to(f"{self.units['energy']} / {self.temperature_unit}").magnitude
    # set radiative wavelength grids
    def _set_wavelength_grids(self):
        self.UV_wavelength_grid = Q_(np.linspace(100.0, 1000.0, 1000), "nanometer")
        self.IR_wavelength_grid = Q_(np.linspace(1000.0, 30000.0, 2000), "nanometer")
    # compute stellar intensity
    def _set_stellar_intensity(self):
        self.stellar_B_lambda = self.stellar_data.B_lambda(self.UV_wavelength_grid)
        plot_stellar_B_lambda(
            wavelength_grid=self.UV_wavelength_grid,
            B_lambda=self.stellar_B_lambda,
            output_file=os.path.join(self.output_dir, "stellar_B_lambda.png"),
        )
    # stellar flux
    def _set_base_stellar_flux(self):
        planet_disk_solid_angle = Q_(np.pi, "steradian")
        self.F0 = (planet_disk_solid_angle * self.stellar_B_lambda *
            (self.stellar_data.radius / self.planet_data.orbital_distance) ** 2
        ).to("W / AU**2 / m")
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
        plot_chemical_mass_profiles(
            mu=mu,
            chem_data=chem_data,
            altitude=z,
            mass_unit=self.mass_unit,
            length_unit="km",
            output_file=os.path.join(self.output_dir, "chemical_mass_profiles_initial_guess.png"),
        )
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
        z_max = self.atmosphere_data["z_max"].to(self.length_unit).magnitude
        # Internal solver variables are unitless floats in the units set by
        # _set_internal_units().
        return np.linspace(
            z_min,
            z_max,
            n_layers,
        )
    # set initial temperature profile
    def _initial_temperature_profile(self, z):
        temperature_surface_guess_K = 288.0
        return np.full_like(z, temperature_surface_guess_K, dtype=float)
    # set initial density
    def _initial_density_profile(self, z):
        density_surface_guess = 1.225    # kg/m^3
        scale_height_guess = 8500.0      # m
        # density profile
        density = density_surface_guess*np.exp(-z/scale_height_guess)
        # compute total mass
        shell_area = 4.0*np.pi*(self.planet_data.planet_radius.to(self.length_unit).magnitude + z)**2
        density_mass = np.trapezoid(density*shell_area, z)
        # renormalize mass
        atmosphere_mass = self._atmosphere_mass.to(self.mass_unit).magnitude
        return density*atmosphere_mass/density_mass
    # gravity profile
    def _gravity_profile(self, z):
        planet_mass = self.planet_data.planet_mass.to(self.mass_unit).magnitude
        planet_radius = self.planet_data.planet_radius.to(self.length_unit).magnitude
        grav_const = G.to(f"{self.length_unit}^3 / {self.mass_unit} / {self.time_unit}^2").magnitude
        return grav_const*planet_mass/(planet_radius + z)**2
    # mean molecular weight profile
    def _mean_molecular_weight_profile(self, chem_data, z, tol=1.e-5):
        n_layers = len(z)
        # avg. molecular weight
        mu_ofz = Q_(np.zeros(n_layers, dtype=float), "amu")
        moles_fraction_sum = np.zeros_like(mu_ofz)
        # run over species / mole_fraction array
        for species, mole_fraction in chem_data.mole_fraction_profiles.items():
            molecular_weight = self._molecular_weight(species)
            xi_ofz = np.asarray(mole_fraction, dtype=float)
            #print(xi_ofz, molecular_weight)
            mu_ofz += xi_ofz*molecular_weight
            moles_fraction_sum += xi_ofz
        if np.any(moles_fraction_sum < 1.0-tol):
            log.warning("Sum_i x_i(z) < 1")
        return mu_ofz.to(self.mass_unit).magnitude
    # molecular weight
    def _molecular_weight(self, species):
        try:
            return Q_(formula(species).mass, "amu")
        except Exception as exc:
            log.error(f"cannot compute molecular weight for species {species}: {exc}")
    # species number density
    def _species_number_density(self, density, chem_data, mu):
        # total number density
        total_number_density = density / mu
        # number density of species
        species_number_density = {}
        for species, mole_fraction in chem_data.mole_fraction_profiles.items():
            xi_ofz = np.asarray(mole_fraction, dtype=float)
            species_number_density[species] = xi_ofz*total_number_density
        return species_number_density
    # surface pressure
    def _set_surface_pressure(self):
        _g = G * self.planet_data.planet_mass / self.planet_data.planet_radius ** 2
        self._P_surf = _g * self._atmosphere_mass / (4. * np.pi * self.planet_data.planet_radius ** 2)
        log.info(f"\t surface pressure: {self._P_surf.to('bar')}")
    # solve hydrostatic pressure profile
    def _solve_hydrostatic_pressure(self, z, temperature, gravity, density=None, mu=None):
        z = np.asarray(z, dtype=float)
        temperature = np.asarray(temperature, dtype=float)
        gravity = np.asarray(gravity, dtype=float)
        # if chemical potential provided
        if mu is not None:
            mu = np.asarray(mu, dtype=float)
            dlogp_dz_abs = mu*gravity/(self._kB*temperature)
            logp_integral = cumulative_trapezoid(dlogp_dz_abs, z, initial=0.0)
            return self._P_surf.to(self.pressure_unit).magnitude*np.exp(-logp_integral)
        # else use density
        if density is None:
            log.error("density profile is required when mu is not provided")
        density = np.asarray(density, dtype=float)
        column_weight = cumulative_trapezoid(density*gravity, z, initial=0.0)
        return self._P_surf.to(self.pressure_unit).magnitude - column_weight
    # compute optical depth
    def _compute_optical_depth(self, config: AtmLayerDyn):
        tau_UV = self._compute_optical_depthUV(config)
        tau_IR = self._compute_optical_depthIR(config)
        return OpticalDepthProfiles(
            tau_UV=tau_UV,
            tau_IR=tau_IR
        )
    # compute UV optical depth
    def _compute_optical_depthUV(self, config: AtmLayerDyn):
        n_wavelength = self.UV_wavelength_grid.magnitude.size
        n_layers = config.altitude.size
        return np.zeros((n_wavelength, n_layers), dtype=float)
    # compute IR optical depth
    def _compute_optical_depthIR(self, config: AtmLayerDyn):
        absorption_coefficient = set_IR_absorption_coefficient_profile(
            wavelength_grid=self.IR_wavelength_grid,
            altitude=config.altitude,
            species_number_density=config.species_number_density,
            temperature=Q_(config.temperature, self.temperature_unit),
            pressure=Q_(config.pressure, self.pressure_unit),
            output_dir=self.output_dir,
        )
        optical_depth = absorption_optical_depth_profile(
            z=config.altitude,
            absorption_coefficient=absorption_coefficient,
        )
        return optical_depth
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
        if self._max_iter_loop < 1:
            log.error("max_iter_loop must be >= 1")
        if self._rel_tol < 0.0 or self._abs_tol < 0.0:
            log.error("loop tolerances must be non-negative")
        if self._damping_loop <= 0.0 or self._damping_loop > 1.0:
            log.error("damping_loop must satisfy 0 < damping_loop <= 1")
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
