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
    def __init__(self, protocell_list):
        self.protocell_list = protocell_list
        # global conditions
        self.temperature = None
        self.pressure = None
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
    # abstract methods
    # ----------------------------
    @abstractmethod
    def update_external_conditions(self, dt):
        raise NotImplementedError("update external conditions")
    @abstractmethod
    def apply_transport(self, dt):
        raise NotImplementedError("apply transport")
    @abstractmethod
    def apply_environmental_effects(self, dt):
        raise NotImplementedError("apply environmaental effects")

