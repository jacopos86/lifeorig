from dataclasses import dataclass


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
