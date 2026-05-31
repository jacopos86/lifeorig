import numpy as np
from scipy.constants import Boltzmann
from src.utilities.logging_module import log


# estimate Bond albedo from absorption/scattering optical depths
def _estimate_bond_albedo(self, config):
    if self.uv_wavelength_grid is None or self.F0 is None:
        log.error("stellar spectrum and flux are required before estimating albedo")
    z = np.asarray(config.altitude, dtype=float)
    pressure = np.asarray(config.pressure, dtype=float)
    temperature = np.asarray(config.temperature, dtype=float)
    mu = np.asarray(config.mean_molecular_mass, dtype=float)
    if pressure.shape != z.shape:
        log.error("pressure profile must have the same shape as altitude grid")
    if temperature.shape != z.shape:
        log.error("temperature profile must have the same shape as altitude grid")
    if mu.shape != z.shape:
        log.error("mean molecular mass profile must have the same shape as altitude grid")
    density = pressure*mu/(Boltzmann*temperature)
    absorption_opacity = _absorption_opacity_profile(self, config)
    scattering_opacity = _scattering_opacity_profile(self, config)
    tau_abs = _vertical_optical_depth(
        self,
        opacity=absorption_opacity,
        density=density,
        z=z,
        name="absorption_opacity",
    )
    tau_sca = _vertical_optical_depth(
        self,
        opacity=scattering_opacity,
        density=density,
        z=z,
        name="scattering_opacity",
    )
    spectral_albedo = _spectral_albedo_from_optical_depth(tau_abs=tau_abs, tau_sca=tau_sca)
    bond_albedo = _flux_weighted_bond_albedo(self, spectral_albedo)
    return bond_albedo, spectral_albedo


# placeholder: fill with species/continuum absorption opacity later
def _absorption_opacity_profile(self, config):
    n_wavelength = self.uv_wavelength_grid.magnitude.size
    n_layers = config.altitude.size
    return np.zeros((n_wavelength, n_layers), dtype=float)


# placeholder: fill with Rayleigh/cloud/aerosol scattering opacity later
def _scattering_opacity_profile(self, config):
    n_wavelength = self.uv_wavelength_grid.magnitude.size
    n_layers = config.altitude.size
    return np.zeros((n_wavelength, n_layers), dtype=float)


# integrate opacity over the vertical column: tau(lambda) = int kappa*rho dz
def _vertical_optical_depth(self, opacity, density, z, name):
    opacity = np.asarray(opacity, dtype=float)
    density = np.asarray(density, dtype=float)
    z = np.asarray(z, dtype=float)
    expected_shape = (self.uv_wavelength_grid.magnitude.size, z.size)
    if opacity.shape != expected_shape:
        log.error(f"{name} must have shape (n_wavelength, n_layers)")
    if np.any(opacity < 0.0):
        log.error(f"{name} must be non-negative")
    if np.any(density < 0.0):
        log.error("density must be non-negative for optical depth integration")
    return np.trapezoid(opacity*density[None, :], z, axis=1)


# two-stream-inspired spectral albedo from absorption/scattering optical depths
def _spectral_albedo_from_optical_depth(tau_abs, tau_sca):
    tau_abs = np.asarray(tau_abs, dtype=float)
    tau_sca = np.asarray(tau_sca, dtype=float)
    if tau_abs.shape != tau_sca.shape:
        log.error("absorption and scattering optical depths must have the same shape")
    if np.any(tau_abs < 0.0) or np.any(tau_sca < 0.0):
        log.error("optical depths must be non-negative")
    tau_tot = tau_abs + tau_sca
    single_scattering_albedo = np.zeros_like(tau_tot)
    has_extinction = tau_tot > 0.0
    single_scattering_albedo[has_extinction] = tau_sca[has_extinction]/tau_tot[has_extinction]
    reflectance_inf = (
        (1.0 - np.sqrt(1.0 - single_scattering_albedo))
        /(1.0 + np.sqrt(1.0 - single_scattering_albedo))
    )
    spectral_albedo = reflectance_inf*(1.0 - np.exp(-tau_tot))
    return np.clip(spectral_albedo, 0.0, 1.0)


# flux-weight spectral albedo into a Bond albedo
def _flux_weighted_bond_albedo(self, spectral_albedo):
    spectral_albedo = np.asarray(spectral_albedo, dtype=float)
    wavelength_nm = self.uv_wavelength_grid.to("nanometer").magnitude
    flux = self.F0.magnitude
    if spectral_albedo.shape != wavelength_nm.shape:
        log.error("spectral albedo must have the same shape as wavelength grid")
    if flux.shape != wavelength_nm.shape:
        log.error("stellar flux must have the same shape as wavelength grid")
    total_flux = np.trapezoid(flux, wavelength_nm)
    if total_flux <= 0.0:
        log.error("integrated stellar flux must be > 0 to estimate Bond albedo")
    bond_albedo = np.trapezoid(spectral_albedo*flux, wavelength_nm)/total_flux
    return float(np.clip(bond_albedo, 0.0, 1.0))
