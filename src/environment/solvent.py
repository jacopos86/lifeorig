from dataclasses import dataclass
from src.common.units import Q_
from src.environment.external_drive_params import LiquidLevelParams

@dataclass
class SolventData:
    name: str
    liquid_level_params: LiquidLevelParams
    density: Q_ | None = None
    dynamic_viscosity: Q_ | None = None
    dielectric_constant: float | None = None
    diffusion_scale: float = 1.0
    polarity: float | None = None