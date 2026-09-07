from dataclasses import dataclass
import numpy as np
import os
from src.common.units import Q_
from src.atmosphere_solver.atmosph_data import AtmLayerDyn
from src.utilities.plot_stellar_spectrum import plot_stellar_B_lambda

#
#   radiative solver
#   1) net energy flux
#   2) update T(z)
#

@dataclass
class RadiativeSolverSettings:
    max_iter: int
    abs_tol: float
    damping: float

class RadiativeSolver:
    def __init__(self, atmosphere_data, stellar_data, planet_data, units, output_dir):
        # store initial data
        self.atmosphere_data = atmosphere_data
        self.stellar_data = stellar_data
        self.planet_data = planet_data
        self._units = units
        self.output_dir = output_dir
        # settings
        radiative_data = atmosphere_data.get("radiative_solver", {}) or {}
        self.solver_type = radiative_data.get("type", "two_stream")
        settings_data = radiative_data.get("settings", {}) or {}
        self.settings = RadiativeSolverSettings(
            max_iter=int(settings_data.get("max_iter", 100)),
            abs_tol=float(settings_data.get("abs_tol", 1.0e-3)),
            damping=float(settings_data.get("damping", 0.5)),
        )
        # internal variables
        self.UV_wavelength_grid = None
        self.IR_wavelength_grid = None
        self.stellar_B_lambda = None
        self.F0 = None
        # 1) set wavelength grids
        self._set_wavelength_grids()
        # 2) set stellar flux
        self._set_stellar_intensity()
        # 3) set base stellar flux
        self._set_base_stellar_flux()
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
