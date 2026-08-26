import numpy as np
from dataclasses import dataclass, field
from src.common.units import Q_

#
#  temporal grid 
#

@dataclass
class TimeGrid:
    # internal unit of time : day
    T: object
    dt: object
    nt: int
    start: object = field(default_factory=lambda: Q_(0.0, "day"))
    time: Q_ = field(default_factory=lambda: np.array([]))
    def __post_init__(self):
        self._validate_dimensions()
        self._normalize_units()
        self._validate_values()
        self._validate_size_consistency()
        self._build_time_array()
    def _validate_dimensions(self):
        if not self.T.check("[time]"):
            raise ValueError("T must have dimensions of time")
        if not self.dt.check("[time]"):
            raise ValueError("dt must have dimensions of time")
        if not self.start.check("[time]"):
            raise ValueError("start must have dimensions of time")
        if not isinstance(self.nt, int):
            raise TypeError("nt must be integer")
    def _normalize_units(self):
        self.T = self.T.to("day")
        self.dt = self.dt.to("day")
        self.start = self.start.to("day")
    def _validate_values(self):
        if self.T.magnitude <= 0:
            raise ValueError("T must be positive")
        if self.dt.magnitude <= 0:
            raise ValueError("dt must be positive")
        if self.nt <= 0:
            raise ValueError("nt must be positive")
    def _validate_size_consistency(self):
        exp_T = self.dt.magnitude * self.nt
        if not (0.999 * exp_T <= self.T.magnitude <= 1.001 * exp_T):
            raise ValueError(
                f"Inconsistent TimeGrid: T={self.T.magnitude} day, "
                f"dt={self.dt.magnitude} day, nt={self.nt}, "
                f"expected total T ≈ {exp_T} day"
            )
    # ----------------------------------
    #    time array
    # ----------------------------------
    def _build_time_array(self):
        """ Build time array in day """
        t_vals = np.linspace(
            self.start.magnitude,
            self.start.magnitude + self.T.magnitude,
            self.nt,
            endpoint=False
        )
        self.time = Q_(t_vals, "day")
    # ---------------------------------
    #    interface
    # ---------------------------------
    def get_time(self, units) -> np.ndarray:
        return self.time.to(units).magnitude
    def __getitem__(self, idx):
        return self.time[idx]
    def __len__(self):
        return self.nt