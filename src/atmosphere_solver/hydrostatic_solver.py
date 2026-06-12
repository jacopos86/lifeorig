from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
from typing import Callable
from src.atmosphere_solver.atm_struct_driver import AtmLayerDyn
from src.common.units import Q_
from src.exo_chem.easy_chem_driver import run_easy_chem_full_profile
from src.utilities.logging_module import log

@dataclass
class PressureDensitySolverSettings:
    max_iter: int
    abs_tol: float
    rel_tol: float
    logp_tol: float
    damping: float = 0.5
    min_pressure: float = 1e-12
    anderson_depth: int = 5

@dataclass
class PressureDensitySolveResult:
    layers: AtmLayerDyn
    converged: bool
    iterations: int
    max_abs_pressure_delta: float
    max_rel_pressure_delta: float
    max_log_pressure_delta: float


MeanMolecularWeightCallback = Callable[[object, np.ndarray], np.ndarray]
SpeciesNumberDensityCallback = Callable[[np.ndarray, object, np.ndarray], dict[str, np.ndarray]]
HydrostaticPressureCallback = Callable[..., np.ndarray]


class FixedTemperaturePressureDensitySolver(ABC):
    """
    Base solver for fixed-temperature pressure/density equilibrium.

    This class intentionally knows only about the pressure-density fixed point.
    The parent atmosphere solver should provide callbacks for hydrostatic
    pressure, mean molecular mass, and species number density so this module
    stays separate from radiative transfer and atmosphere setup details.
    """

    solver_name = "fixed-temperature pressure-density"

    def __init__(
        self,
        *,
        settings: PressureDensitySolverSettings,
        pressure_unit: str,
        temperature_unit: str,
        density_unit: str,
        kB: float,
        atomic_abundances: dict[str, float],
        chemical_species,
        solve_hydrostatic_pressure: HydrostaticPressureCallback,
        mean_molecular_weight_profile: MeanMolecularWeightCallback,
        species_number_density: SpeciesNumberDensityCallback,
    ):
        self.settings = settings
        self.pressure_unit = pressure_unit
        self.temperature_unit = temperature_unit
        self.density_unit = density_unit
        self._kB = kB
        self.atomic_abundances = atomic_abundances
        self.chemical_species = chemical_species
        self._solve_hydrostatic_pressure = solve_hydrostatic_pressure
        self._mean_molecular_weight_profile = mean_molecular_weight_profile
        self._species_number_density = species_number_density

    def solve(self, config: AtmLayerDyn) -> PressureDensitySolveResult:
        z, temperature, pressure, gravity = self._unpack_config(config)
        log_pressure = np.log(np.maximum(pressure, self.settings.min_pressure))

        converged = False
        iterations = 0
        metrics = {
            "max_abs": np.inf,
            "max_rel": np.inf,
            "max_dlogp": np.inf,
        }

        self._start()

        for iteration in range(self.settings.max_iter):
            pressure_old = np.exp(log_pressure)
            hydro_pressure, chem_data, mu = self._fixed_point_pressure(
                z=z,
                temperature=temperature,
                gravity=gravity,
                pressure=pressure_old,
            )
            residual = np.log(hydro_pressure) - log_pressure
            metrics = self._pressure_metrics(
                pressure_old=pressure_old,
                pressure_new=hydro_pressure,
                residual=residual,
            )

            iterations = iteration + 1
            self._log_iteration(iterations, metrics)

            if self._has_converged(metrics):
                converged = True
                log_pressure = np.log(hydro_pressure)
                log.info(f"{self.solver_name} converged in {iterations} iterations")
                break

            log_pressure = self._next_log_pressure(
                iteration=iteration,
                log_pressure=log_pressure,
                residual=residual,
            )

        if not converged:
            log.warning(f"{self.solver_name} did not converge after {self.settings.max_iter} iterations")

        pressure = np.maximum(np.exp(log_pressure), self.settings.min_pressure)
        layers = self._final_layers(
            z=z,
            temperature=temperature,
            pressure=pressure,
            gravity=gravity,
        )

        return PressureDensitySolveResult(
            layers=layers,
            converged=converged,
            iterations=iterations,
            max_abs_pressure_delta=metrics["max_abs"],
            max_rel_pressure_delta=metrics["max_rel"],
            max_log_pressure_delta=metrics["max_dlogp"],
        )

    def _unpack_config(self, config: AtmLayerDyn):
        z = np.asarray(config.altitude, dtype=float)
        temperature = np.asarray(config.temperature, dtype=float)
        pressure = np.asarray(config.pressure, dtype=float)
        gravity = np.asarray(config.gravity, dtype=float)

        if np.any(temperature <= 0.0):
            log.error("temperature must be positive for pressure-density equilibrium solve")
        if np.any(pressure <= 0.0):
            log.error("pressure must be positive for pressure-density equilibrium solve")
        return z, temperature, pressure, gravity

    def _start(self):
        log.warning(
            f"Starting {self.solver_name}: maxiter={self.settings.max_iter}, "
            f"abs_tol={self.settings.abs_tol}, rel_tol={self.settings.rel_tol}, "
            f"logp_tol={self.settings.logp_tol}"
        )

    def _fixed_point_pressure(self, *, z, temperature, gravity, pressure):
        chem_data = run_easy_chem_full_profile(
            pressure=Q_(pressure, self.pressure_unit),
            temperature=Q_(temperature, self.temperature_unit),
            atomic_abund=self.atomic_abundances,
            chemical_species=self.chemical_species,
        )
        mu = self._mean_molecular_weight_profile(chem_data, z)
        if np.any(mu <= 0.0):
            log.error("mean molecular mass must be positive for pressure-density equilibrium solve")

        hydro_pressure = self._solve_hydrostatic_pressure(
            z=z,
            temperature=temperature,
            gravity=gravity,
            mu=mu,
        )
        hydro_pressure = np.maximum(hydro_pressure, self.settings.min_pressure)
        return hydro_pressure, chem_data, mu

    def _pressure_metrics(self, *, pressure_old, pressure_new, residual):
        pressure_delta = np.abs(pressure_new - pressure_old)
        pressure_scale = np.maximum(np.maximum(np.abs(pressure_new), np.abs(pressure_old)), 1e-30)
        return {
            "max_abs": float(np.max(pressure_delta)),
            "max_rel": float(np.max(pressure_delta / pressure_scale)),
            "max_dlogp": float(np.max(np.abs(residual))),
        }

    def _has_converged(self, metrics: dict[str, float]) -> bool:
        return (
            metrics["max_dlogp"] <= self.settings.logp_tol
            or metrics["max_abs"] <= self.settings.abs_tol
            or metrics["max_rel"] <= self.settings.rel_tol
        )

    def _log_iteration(self, iteration: int, metrics: dict[str, float]):
        log.warning(
            f"{self.solver_name} iter {iteration}: "
            f"max_abs_dP={metrics['max_abs']:.3e} Pa, "
            f"max_rel_dP={metrics['max_rel']:.3e}, "
            f"max_dlogP={metrics['max_dlogp']:.3e}"
        )

    def _final_layers(self, *, z, temperature, pressure, gravity) -> AtmLayerDyn:
        chem_data = run_easy_chem_full_profile(
            pressure=Q_(pressure, self.pressure_unit),
            temperature=Q_(temperature, self.temperature_unit),
            atomic_abund=self.atomic_abundances,
            chemical_species=self.chemical_species,
        )
        mu = self._mean_molecular_weight_profile(chem_data, z)
        density = pressure * mu / (self._kB * temperature)
        species_number_density = self._species_number_density(
            density=density,
            chem_data=chem_data,
            mu=mu,
        )
        return AtmLayerDyn(
            altitude=z,
            pressure=pressure,
            temperature=temperature,
            mean_molecular_mass=mu,
            gravity=gravity,
            species_number_density=species_number_density,
            chemistry=chem_data,
        )

    @abstractmethod
    def _next_log_pressure(self, *, iteration: int, log_pressure: np.ndarray, residual: np.ndarray) -> np.ndarray:
        raise NotImplementedError


class StandardPressureDensitySolver(FixedTemperaturePressureDensitySolver):
    solver_name = "standard fixed-temperature pressure-density solver"

    def _next_log_pressure(self, *, iteration: int, log_pressure: np.ndarray, residual: np.ndarray) -> np.ndarray:
        damping = max(0.0, min(1.0, self.settings.damping))
        return log_pressure + damping * residual


class AndersonPressureDensitySolver(FixedTemperaturePressureDensitySolver):
    solver_name = "Anderson fixed-temperature pressure-density solver"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._residual_history: list[np.ndarray] = []

    def _start(self):
        self._residual_history.clear()
        super()._start()
        log.info(f"Anderson history depth={self.settings.anderson_depth}")

    def _next_log_pressure(self, *, iteration: int, log_pressure: np.ndarray, residual: np.ndarray) -> np.ndarray:
        self._residual_history.append(residual.copy())

        if len(self._residual_history) == 1:
            return log_pressure + self.settings.damping * residual

        depth = min(self.settings.anderson_depth, len(self._residual_history) - 1)
        residual_basis = np.column_stack(self._residual_history[-depth:])

        try:
            alpha, _, _, _ = np.linalg.lstsq(residual_basis, residual, rcond=None)
        except np.linalg.LinAlgError:
            log.warning("Anderson acceleration failed; falling back to damped update")
            return log_pressure + self.settings.damping * residual

        correction = residual_basis @ alpha
        log.info(f"Anderson acceleration applied with depth={depth}")
        return log_pressure + correction
