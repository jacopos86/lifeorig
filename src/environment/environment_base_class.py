from abc import ABC, abstractmethod

#
#   Environment base class
#

class Environment(ABC):
    """
    Abstract environment class.

    Handles:
    - external conditions (T, gradients, chemistry)
    - protocell placement
    - transport between protocells and environment
    """
    def __init__(self):
        # global conditions
        self.temperature = None
        self.pressure = None
        self.g = None
        # external molecule reservoir (optional)
        self.external_field = {}
    # ----------------------------
    # main evolution step
    # ----------------------------
    def step(self, dt):
        """
        Advance environment + protocell coupling
        """
        self.update_external_conditions(dt)
        self.apply_transport(dt)
        self.apply_environmental_effects(dt)
    # ----------------------------
    # get environment state
    # ----------------------------
    def get_environment_state(self):
        state = self._build_state()
        return {
            "time": self.state.time,
            "liquid_level": self.state.liquid_level,
            "temperature": self.state.temperature,
            "pressure": self.state.pressure,
            "gravity": self.state.gravity,
            "pH": self.state.pH,
            "ionic_strength": self.state.ionic_strength,
            "uv_flux": self.state.uv_flux,
            "sediment_surface_area": self.state.sediment_surface_area,
            "surface_catalytic_probability_z": self.surface_catalytic_probability_z,
            "volume": self.state.volume
        }
    # ----------------------------
    # abstract methods
    # ----------------------------
    @abstractmethod
    def update_external_conditions(self, dt):
        raise NotImplementedError("update external conditions")
    @abstractmethod
    def set_spatial_geometry(self, geometry):
        raise NotImplementedError("set the spatial geometry")