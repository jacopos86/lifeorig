from abc import ABC, abstractmethod
from dataclasses import dataclass
from src.common.units import Q_

# =========================================================
# parameter containers
# =========================================================

@dataclass
class LiquidLevelParams:
    model_type: str                  # "constant", "sinusoidal", "piecewise"
    base_level: Q_
    amplitude: Q_ | None
    period: Q_ | None
    phase: float = 0.0
    switch_times: Q_ | None = None
    levels: Q_ | None = None


@dataclass
class RadiationParams:
    model_type: str                  # "constant", "sinusoidal", "pulse"
    base_level: float
    amplitude: float = 0.0
    period: float = 1.0
    phase: float = 0.0
    pulse_start: float | None = None
    pulse_end: float | None = None
    pulse_level: float | None = None


# =========================================================
# abstract base classes
# =========================================================

class TimeDependentField(ABC):
    def __init__(self, params):
        self.params = params
    @abstractmethod
    def value(self, t: Q_) -> Q_:
        raise NotImplementedError
