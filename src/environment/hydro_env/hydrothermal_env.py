from src.environment.environment_base_class import Environment

#
#    Hydrothermal vent (1D) for simplicity
#

class HydrothermalVent(Environment):
    def __init__(self, protocell_list):
        super().__init__(protocell_list)
        # define gradient
        self.T_hot = None    # near vent
        self.T_cold = None   # far field
    def update_external_conditions(self, dt):
        # temperature depends on position
        for i, proto in enumerate(self.protocell_list):
            x = self.positions[i]
            proto.temperature = self.T_hot * (1 - x) + self.T_cold * x
    def apply_transport(self, dt):
        # simple exchange with environment
        for proto in self.protocell_list:
            for mol, metab in proto.M_set.items():
                influx = 0.01 * dt  # toy model
                metab.count += influx
    def apply_environmental_effects(self, dt):
        # temperature-dependent reaction scaling
        for proto in self.protocell_list:
            T = proto.temperature
            proto.reaction_scale = np.exp(-1.0 / T)  # placeholder