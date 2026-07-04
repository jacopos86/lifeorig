from dataclasses import dataclass, field
import numpy as np
from src.environment.environment_base_class import Environment

# 
#   Vent cell
#

@dataclass
class VentCell:
    """
    One coarse slice of the porous hydrothermal chimney wall.

    The cell is not empty ocean water. It represents water-filled pores plus
    mineral wall surface in a small volume of vent rock.
    """
    idx: int
    x: float
    temperature: float
    pressure: float
    pH: float
    flow_velocity: float
    porosity: float
    bulk_volume: float
    fluid_volume: float
    mineral_surface_area: dict[str, float] = field(default_factory=dict)
    concentrations: dict[str, float] = field(default_factory=dict)

#
#   Pore compartment
#

@dataclass
class PoreCompartment:
    """
    Optional explicit pore reactor inside one VentCell.

    First simulations can ignore these and use only VentCell-averaged porosity
    and mineral surface area.
    """
    idx: int
    cell_idx: int
    radius: float
    length: float
    mineral_type: str
    concentrations: dict[str, float] = field(default_factory=dict)
    @property
    def volume(self):
        return np.pi * self.radius**2 * self.length
    @property
    def surface_area(self):
        return 2.0 * np.pi * self.radius * self.length
    @property
    def residence_time(self):
        return None

#
#    Hydrothermal vent: 1D porous mineral reactor
#

class HydrothermalVent(Environment):
    """
    Coarse-grained 1D hydrothermal vent environment.

    Geometry:
        hot vent fluid -> porous mineral chimney wall -> cold ocean water

    The 1D grid is a stack of porous rock cells. Each cell contains pore water,
    dissolved species, mineral surface area, flow, temperature, and pH.
    """
    def __init__(self, protocell_list=None, data=None):
        super().__init__(protocell_list or [])
        data = data or {}
        self.n_cells = int(data.get("n_cells", 10))
        self.length = self._as_float(data.get("length"), default=1.0)
        self.cross_section_area = self._as_float(data.get("cross_section_area"), default=1.0)
        self.T_hot = self._as_float(data.get("T_hot", data.get("temperature_hot")), default=400.0)
        self.T_cold = self._as_float(data.get("T_cold", data.get("temperature_cold")), default=280.0)
        self.pressure = self._as_float(data.get("pressure"), default=1.0e7)
        self.pH_hot = float(data.get("pH_hot", 10.5))
        self.pH_cold = float(data.get("pH_cold", 8.0))
        self.flow_velocity = self._as_float(data.get("flow_velocity"), default=1.0e-5)
        self.porosity = float(data.get("porosity", 0.3))
        self.surface_area_per_volume = self._as_float(
            data.get("surface_area_per_volume"),
            default=1.0e5,
        )
        self.vent_concentrations = dict(data.get("vent_concentrations", {}))
        self.ocean_concentrations = dict(data.get("ocean_concentrations", {}))
        self.mineral_fractions = dict(
            data.get(
                "mineral_fractions",
                {
                    "FeS": 0.5,
                    "NiS": 0.1,
                    "magnetite": 0.2,
                    "brucite": 0.1,
                    "silica": 0.1,
                },
            )
        )
        self.cells = self._build_cells()
        self.pores = self._build_pores(data)
        self.coacervates = []
        self.vesicle_candidates = []
        self.time = 0.0
    # build internal cells
    def _build_cells(self):
        cells = []
        dx = self.length / self.n_cells
        bulk_volume = self.cross_section_area * dx
        fluid_volume = self.porosity * bulk_volume
        total_surface_area = self.surface_area_per_volume * bulk_volume

        for idx in range(self.n_cells):
            if self.n_cells == 1:
                xfrac = 0.0
            else:
                xfrac = idx / (self.n_cells - 1)

            mineral_surface_area = {
                mineral: frac * total_surface_area
                for mineral, frac in self.mineral_fractions.items()
            }
            cells.append(
                VentCell(
                    idx=idx,
                    x=xfrac * self.length,
                    temperature=self._linear(self.T_hot, self.T_cold, xfrac),
                    pressure=self.pressure,
                    pH=self._linear(self.pH_hot, self.pH_cold, xfrac),
                    flow_velocity=self.flow_velocity,
                    porosity=self.porosity,
                    bulk_volume=bulk_volume,
                    fluid_volume=fluid_volume,
                    mineral_surface_area=mineral_surface_area,
                    concentrations=self._mixed_concentrations(xfrac),
                )
            )
        return cells
    # build pores
    def _build_pores(self, data):
        pores = []
        n_pores_per_cell = int(data.get("n_pores_per_cell", 0))
        if n_pores_per_cell <= 0:
            return pores

        radius = self._as_float(data.get("pore_radius"), default=1.0e-6)
        length = self._as_float(data.get("pore_length"), default=self.length / self.n_cells)
        minerals = list(self.mineral_fractions) or ["FeS"]
        pore_idx = 0

        for cell in self.cells:
            for ipore in range(n_pores_per_cell):
                mineral_type = minerals[ipore % len(minerals)]
                pores.append(
                    PoreCompartment(
                        idx=pore_idx,
                        cell_idx=cell.idx,
                        radius=radius,
                        length=length,
                        mineral_type=mineral_type,
                        concentrations=dict(cell.concentrations),
                    )
                )
                pore_idx += 1
        return pores

    def _mixed_concentrations(self, xfrac):
        species = set(self.vent_concentrations) | set(self.ocean_concentrations)
        return {
            species_name: self._linear(
                float(self.vent_concentrations.get(species_name, 0.0)),
                float(self.ocean_concentrations.get(species_name, 0.0)),
                xfrac,
            )
            for species_name in species
        }

    def step(self, dt):
        self.time += dt
        self.update_external_conditions(dt)
        self.apply_transport(dt)
        self.apply_environmental_effects(dt)

    def update_external_conditions(self, dt):
        return

    def apply_transport(self, dt):
        return

    def apply_environmental_effects(self, dt):
        return

    def get_cell(self, cell_idx):
        if len(self.cells) == 0:
            raise RuntimeError("HydrothermalVent has no cells")
        return self.cells[cell_idx]

    def get_cell_at_position(self, x):
        if len(self.cells) == 0:
            raise RuntimeError("HydrothermalVent has no cells")
        x = max(0.0, min(self.length, x))
        xfrac = x / self.length if self.length > 0.0 else 0.0
        idx = int(round(xfrac * (len(self.cells) - 1)))
        idx = max(0, min(len(self.cells) - 1, idx))
        return self.cells[idx]

    def get_environment_state(self, cell_idx):
        cell = self.get_cell(cell_idx)
        return {
            "temperature": cell.temperature,
            "pressure": cell.pressure,
            "pH": cell.pH,
            "flow_velocity": cell.flow_velocity,
            "porosity": cell.porosity,
            "fluid_volume": cell.fluid_volume,
            "mineral_surface_area": cell.mineral_surface_area,
            "concentrations": cell.concentrations,
        }

    @staticmethod
    def _linear(left, right, xfrac):
        return (1.0 - xfrac) * left + xfrac * right

    @staticmethod
    def _temperature_scale(temperature):
        return float(np.exp(-1.0 / max(temperature, 1.0)))

    @staticmethod
    def _as_float(value, default):
        if value is None:
            return float(default)
        if hasattr(value, "to"):
            try:
                return float(value.to_base_units().magnitude)
            except Exception:
                return float(value.magnitude)
        return float(value)