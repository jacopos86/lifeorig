from abc import ABC, abstractmethod
from dataclasses import dataclass
import math
from src.common.units import Q_

@dataclass
class ExternalDriveParams:
    model_type: str
    base_level: Q_
    amplitude: Q_ | None = None
    period: Q_ | None = None
    phase: float = 0.0
    switch_times: Q_ | None = None
    levels: Q_ | None = None

# =========================================================
# time profiles for external forces
# =========================================================

class TimeDependentField(ABC):
    def __init__(self, params):
        self.params = params
    @abstractmethod
    def value(self, t: Q_) -> Q_:
        raise NotImplementedError

class Constant(TimeDependentField):
    def value(self, t: Q_) -> Q_:
        return self.params.base_level

class Sinusoidal(TimeDependentField):
    def value(self, t: Q_) -> Q_:
        phase = (2.0 * math.pi * t / self.params.period).to_base_units().magnitude + self.params.phase
        return (
            self.params.base_level
            + self.params.amplitude
            * math.sin(phase)
        )

class PieceWise(TimeDependentField):
    def value(self, t: Q_) -> Q_:
        for idx, switch_time in enumerate(self.params.switch_times):
            if t < switch_time:
                return self.params.levels[idx]
        return self.params.levels[-1]

# =========================================================
# build external drive
# =========================================================

def build_external_drive(params: ExternalDriveParams):
    model_type = params.model_type.lower()
    if model_type == "constant":
        return Constant(params)
    if model_type == "sinusoidal":
        return Sinusoidal(params)
    if model_type == "piecewise":
        return PieceWise(params)
    raise ValueError(f"Unknown external drive model_type: {params.model_type}")
