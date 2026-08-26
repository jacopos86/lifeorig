from dataclasses import dataclass, field
import numpy as np
from src.environment.environment_base_class import Environment
from src.grids.spatial_grid import uniform_spatial_grid_1d
from src.utilities.plot_titan_pool_profile import titan_surface_area_profile
from src.environment.liquid_level_model import Liquid
from src.common.units import Q_

@dataclass
class PondState:
    """
    Well-mixed shallow liquid pond state.

    The liquid is generic: it can be water, brine, methane/ethane, or another
    configured solvent. Aqueous quantities such as pH and ionic strength are
    optional because they are not meaningful for every solvent.
    """
    time: float
    liquid_level: float
    temperature: float
    pressure: float
    volume: float
    gravity: object | None = None
    pH: float | None = None
    ionic_strength: float | None = None
    uv_flux: float | None = None
    sediment_surface_area: float | None = None
    concentrations: dict[str, float] = field(default_factory=dict)


class SurfacePondEnvironment(Environment):
    """
    Generic surface pond environment.

    Geometry:
        one well-mixed shallow liquid body above a sediment/mineral or reactive
        surface

    Main drivers:
        liquid-level cycling, radiation exposure, atmospheric input, reactive
        surface chemistry, concentration during evaporation, dilution during
        refilling, gravity, thermal forcing, and later hydrodynamic mixing.

    The object should consume planet, star, and parsed local environment data,
    then derive the internal parameters needed by the simulation. It should not
    assume an aqueous solvent; pH is optional and should only be used for
    solvent systems where it is physically defined.
    """
    def __init__(
        self,
        planet=None,
        star=None,
        local_env_data=None,
        time_grid=None,
        working_dir=None
    ):
        super().__init__()
        # internal parameters
        self.planet = planet
        self.star = star
        data = local_env_data or {}
        self.local_env_data = data
        self.time = 0.0
        # pond geometry
        geometry = data.get("pool_geometry", {})
        self.set_spatial_geometry(geometry)
        # gravity
        self.g = self.planet.surface_gravity
        # liquid state 
        solvent_data = data.get("solvent_data")
        self.liquid = Liquid(
            solvent_data=solvent_data,
            liquid_level_params=data.get("liquid_level_params"),
            external_forces=data.get("external_forces", {}),
            atmospheric_composition=self.planet.chemical_stationary_config.atmosph_composition(z=self.grid.cell_centers[0])
        )
        # plot liquid state over time
        self.liquid.plot_liquid_level(time_grid=time_grid, working_dir=working_dir)
        # set UV radiation flux
        exit()
        self.temperature = data.get("temperature")
        self.pressure = data.get("pressure")
        self.pH = None if data.get("pH") is None else float(data.get("pH"))
        self.ionic_strength = None if data.get("ionic_strength") is None else float(data.get("ionic_strength"))
        self.uv_flux = float(data.get("uv_flux", 1.0))
        self.atmospheric_input = dict(data.get("atmospheric_input", {}))
        self.washout_rate = float(data.get("washout_rate", 0.0))
        self.photolysis_rates = dict(data.get("photolysis_rates", {}))

        self.state = self._build_state()

    def update_external_conditions(self, dt):
        self.liquid.step(dt=dt, time=self.time, surface_area=self.area)
        self.time += dt
        self.state = self._build_state()

    def apply_transport(self, dt):
        self._apply_sources_and_losses(dt)

    def apply_environmental_effects(self, dt):
        return

    def _build_state(self):
        return PondState(
            time=self.time,
            liquid_level=self.liquid.level,
            temperature=self.temperature,
            pressure=self.pressure,
            volume=self.area * self.liquid.level,
            gravity=self.g,
            pH=self.pH,
            ionic_strength=self.ionic_strength,
            uv_flux=self.uv_flux,
            sediment_surface_area=self.sediment_surface_area,
            concentrations=dict(self.liquid.concentrations),
        )

    def _apply_sources_and_losses(self, dt):
        for species, source_rate in self.atmospheric_input.items():
            self.liquid.concentrations[species] = self.liquid.concentrations.get(species, 0.0) + source_rate * dt
        for species, rate in self.photolysis_rates.items():
            if species in self.liquid.concentrations:
                self.liquid.concentrations[species] *= np.exp(-rate * self.uv_flux * dt)
        if self.washout_rate > 0.0:
            loss = np.exp(-self.washout_rate * dt)
            for species in self.liquid.concentrations:
                self.liquid.concentrations[species] *= loss

    # set grid geometry + surface catalytic activity
    def set_spatial_geometry(self, geometry):
        # set geometry
        self.geometry = geometry
        height = geometry.get("height")
        surface_area_z0 = geometry.get("surface_area_z0")
        n_cells = int(geometry.get("n_grid_cells", 100))
        if height is None:
            raise ValueError("Missing pool_geometry.height")
        if surface_area_z0 is None:
            raise ValueError("Missing pool_geometry.surface_area_z0")
        # vertical 1D grid only
        self.grid = uniform_spatial_grid_1d(length=height, n_cells=n_cells)
        # set catalytic surface probability
        self.reactive_area_to_surface_area_ratio = float(
            geometry.get("reactive_area_to_surface_area_ratio", 0.0)
        )
        probability_hitting_surface_z = self._surface_hit_probability_z(
            surface_area_profile=geometry.get("surface_area_profile"),
            surface_area_z0=surface_area_z0,
            height=height,
        )
        self.surface_catalytic_probability_z = np.clip(
            self.reactive_area_to_surface_area_ratio * probability_hitting_surface_z,
            0.0,
            1.0,
        )
    # set surface hit probability(z)
    def _surface_hit_probability_z(self, surface_area_profile, surface_area_z0, height):
        surface_area_z = titan_surface_area_profile(
            z=self.grid.cell_centers,
            height=height,
            surface_area=surface_area_z0,
            profile_data=surface_area_profile,
        )
        return (surface_area_z / surface_area_z0).to_base_units().magnitude
