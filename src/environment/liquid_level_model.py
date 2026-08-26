from __future__ import annotations
from dataclasses import dataclass
import os
import numpy as np
import matplotlib.pyplot as plt
from src.common.units import Q_
from src.environment.external_drive_params import (
    Constant,
    PieceWise,
    Sinusoidal,
    TimeDependentField,
    build_external_drive,
)
from src.utilities.logging_module import log
from src.environment.solvent import Solvent

#
#   Liquid level parameters
#

@dataclass
class LiquidLevelParams:
    base_level: Q_
    model_type: str = "constant"
    min_level: Q_ | None = None
    max_level: Q_ | None = None
    initial_level: Q_ | None = None

#
#   build liquid level field
#

def build_liquid_level_field(params: LiquidLevelParams) -> TimeDependentField:
    model_type = params.model_type.lower()
    if model_type == "constant":
        return Constant(params)
    if model_type == "sinusoidal":
        return Sinusoidal(params)
    if model_type == "piecewise":
        return PieceWise(params)
    raise log.error(f"Unknown liquid level model_type: {params.model_type}")

#
#   derive chemical water factor
#

def derive_liquid_base_factor(
        env_model: str,
        env_data: dict,
        planet_data: PlanetaryEnvironmentParams,
        chemical_env: ChemEnvResult,
        vesc: Q_,
        min_factor: float = 0.0,
        max_factor: float = 2.0
    ) -> float:
    # check env model
    if env_model == "volcanic_rock":
        return derive_volcanic_rock_water_base_factor(
            chemical_env=chemical_env,
            vesc=vesc,
            min_factor=min_factor,
            max_factor=max_factor
        )
    if env_model == "hydro_vent":
        return None
    if env_model == "Titan":
        return derive_Titan_solvent_base_factor(
            planet_data=planet_data,
            env_data=env_data
        )
    log.error(f"Unknown environment model: {env_model}")

#
#    liquid amplitude factor
#

def derive_liquid_amplitude_factor(
        env_model: str,
        planet_data: PlanetaryEnvironmentParams,
        water_base_factor: float,
        min_factor: float = 0.0,
        max_factor: float = 2.0
    ) -> float:
    if env_model == "hydro_vent":
        return None
    if env_model == "volcanic_rock":
        return derive_volcanic_rock_water_amplitude_factor(
            planet_data=planet_data,
            water_base_factor=water_base_factor,
            min_factor=min_factor,
            max_factor=max_factor
        )
    log.error(f"Unknown environment model: {env_model}")

#
#  Liquid class
#

class Liquid:
    def __init__(
        self,
        solvent_data,
        liquid_level_params=None,
        external_forces=None,
        atmospheric_composition=None
    ):
        # solvent
        self.solvent = Solvent(solvent_data)
        # atmospheric composition
        self.atmospheric_composition = atmospheric_composition or {}
        # set variable liquid levels model
        self._set_liquid_level_model(liquid_level_params, external_forces)
    # set internal liquid model
    def _set_liquid_level_model(self, liquid_level_params, external_forces):
        external_forces = external_forces or {}
        self.level_params = liquid_level_params
        # set initial level
        self.level = (
            liquid_level_params.initial_level
            if liquid_level_params.initial_level is not None
            else liquid_level_params.base_level
        )
        self.rainfall = self._build_external_force(external_forces.get("rainfall"))
        self.evaporation = self._build_external_force(external_forces.get("evaporation"))
        # set initial composition
        self.composition = self.solvent.normalized_composition()
    # build external fields
    def _build_external_force(self, params):
        if params is None:
            return None
        return build_external_drive(params)
    # dynamical step
    def step(self, dt, time, surface_area):
        old_level = self.level
        rain_level = self._level_change(self.rainfall, dt, time)
        evap_level = self._level_change(self.evaporation, dt, time)
        new_level = old_level + rain_level - evap_level
        if self.level_params.min_level is not None and new_level < self.level_params.min_level:
            new_level = self.level_params.min_level
        if self.level_params.max_level is not None and new_level > self.level_params.max_level:
            new_level = self.level_params.max_level
        self._update_concentrations(old_level, new_level, rain_level, evap_level, surface_area)
        self.level = new_level
    # change liquid level
    def _level_change(self, field, dt, time):
        if field is None:
            return 0.0 * self.level.units
        return field.value(time) * dt
    # update concentrations
    def _update_concentrations(self, old_level, new_level, rain_level, evap_level, surface_area):
        if new_level <= 0.0 * new_level.units:
            self.concentrations = {}
            self.composition = self.concentrations
            return
        species_set = set(self.concentrations) | set(self.atmospheric_composition) | self.solvent.species()
        new_concentrations = {}
        for species in species_set:
            old_amount = self.concentrations.get(species, 0.0) * old_level
            evap_amount = self.solvent.mole_fraction(species, 0.0) * evap_level
            rain_amount = self.atmospheric_composition.get(species, 0.0) * rain_level
            new_concentrations[species] = max(((old_amount - evap_amount + rain_amount) / new_level).magnitude, 0.0)
        self.concentrations = new_concentrations
    #
    # INTERFACE
    # plot liquid level
    def plot_liquid_level(self, time_grid, working_dir=None, output_file=None):
        output_file = output_file or os.path.join(working_dir or ".", "liquid_level.pdf")
        level = self.level
        x = time_grid.get_time("day")
        y = []
        # run over time grid
        for time in time_grid.time:
            y.append(level.to(self.level.units).magnitude)
            rain_level = self._level_change(self.rainfall, time_grid.dt, time)
            evap_level = self._level_change(self.evaporation, time_grid.dt, time)
            level = level + rain_level - evap_level
            if self.level_params.min_level is not None and level < self.level_params.min_level:
                level = self.level_params.min_level
            if self.level_params.max_level is not None and level > self.level_params.max_level:
                level = self.level_params.max_level
        fig, ax = plt.subplots()
        ax.plot(x, y)
        ax.set_xlabel("time [day]")
        ax.set_ylabel(f"liquid level [{self.level.units:~P}]")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(output_file)
        plt.close(fig)