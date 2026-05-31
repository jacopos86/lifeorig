from src.environment.environment_base_class import Environment
from src.environment.liquid_level_model import build_liquid_level_field
from src.utilities.logging_module import log
from src.input_data.read_input import p

#
#   Single rock pore
#

class SingleInterstice:
    def __init__(self, idx, height, radius, neighbors=None):
        self.idx = idx
        # structural parameters
        self.height = height
        self.radius = radius
        # neighbors
        self.neighbors = neighbors
        # list of protocells
        self.protocell_list = []
    def is_flooded(self, water_level):
        return water_level >= self.height
    def update_position(self, dt, rng):
        dz = -self.v_settle * dt + np.sqrt(2.0 * self.D_move * dt) * rng.normal()
        self.z += dz
        if self.z <= 0.0:
            self.z = 0.0
            self.surface_contact = True
        else:
            self.surface_contact = False
        if self.z >= pore.height:
            self.z = pore.height
    def info(self):
        log.info(f"\t PORE RADIUS: {self.radius}")
        log.info(f"\t PORE HEIGHT: {self.height}")
        log.info(f"\t PORE NEIGHBORS: {self.neighbors}")
        log.info("\t " + p.sep)
        log.info("\n")

#
#     Porous volcanic rock
#

class VolcanicRock(Environment):
    def __init__(self, protocell_list, data):
        super().__init__(protocell_list)
        # min. surface distance
        self.min_surface_distance = None
        # define mineral catalytic strength
        self.k_surface = None  # set later
        # set list of interstices
        self.interstices = []
        n = data.get("num_pores")
        for i in range(n):
            idx = i+1
            neighbors = []
            # left neighbor
            if i > 0:
                neighbors.append(idx-1)        # i-1 in 0-based → idx-1 in 1-based
            # right neighbor
            if i < n - 1:
                neighbors.append(idx+1)        # i+1 in 0-based → idx+1 in 1-based
            self.interstices.append(
                SingleInterstice(
                    idx=idx,
                    height=data.get("pore_height"),
                    radius=data.get("pore_radius"),
                    neighbors=neighbors
                )
            )
        self.interstices[0].info()
        self.liquid_level = build_liquid_level_field(
            data.get("solvent_data").liquid_level_params
        )
        self.radiation_level = None
        exit()
        # transport parameters
        self.D_solvent = 1.0e-9
        self.gravity_drift = None
        # assign protocells initially
        self._place_initial_protocells()
    def update_external_conditions(self, dt):
        # mostly constant environment
        for proto in self.protocell_list:
            proto.temperature = self.temperature
    def apply_transport(self, dt):
        # weak exchange with environment
        for proto in self.protocell_list:
            for mol, metab in proto.M_set.items():
                metab.count += 0.001 * dt
    def apply_environmental_effects(self, dt):
        # modify outer shell catalysis
        for proto in self.protocell_list:
            if not getattr(proto, "attached_to_surface", False):
                continue
            outer_shell = proto.n_shells - 1
            # attach mineral catalytic effect
            proto.surface_catalysis_active = True
            proto.surface_shell = outer_shell
    def _place_initial_protocells(self):
        if len(self.interstices) == 0:
            log.error("No pores defined")
        exit()
        # randomly assign attachment
        for proto in self.protocell_list:
            proto.attached_to_surface = (np.random.rand() < 0.5)
        for i, proto in enumerate(self.protocell_list):
            j = i % len(self.interstices)
            self.interstices[j].protocells.append(proto)
            proto.interstice_id = self.interstices[j].idx
    def update_water_level(self, dt):
        # placeholder: user-defined forcing
        self.time += dt
        self.water_level = 0.5 * (1.0 + np.sin(self.time))
    def get_flooded_interstices(self):
        return [it for it in self.interstices if it.is_flooded(self.water_level)]
    def step(self, dt):
        self.update_water_level(dt)
        self.mix_connected_interstices(dt)
    def mix_connected_pores(self, dt):
        flooded = self.get_flooded_interstices()
        if len(flooded) < 2:
            return
        moves = []
        for i, src in enumerate(flooded):
            for j, dst in enumerate(flooded):
                if i == j:
                    continue
                rij = np.linalg.norm(src.position - dst.position)
                if rij == 0.0:
                    continue
                # diffusion contribution
                k_diff = self.D_water / (rij * rij)
                # gravity / downhill drift
                dz = src.height - dst.height
                k_grav = self.gravity_drift * max(0.0, dz)
                k_move = k_diff + k_grav
                p_move = min(1.0, k_move * dt)
                for proto in src.protocells:
                    if np.random.rand() < p_move:
                        moves.append((proto, src, dst))
        for proto, src, dst in moves:
            if proto in src.protocells:
                src.protocells.remove(proto)
                dst.protocells.append(proto)
                proto.interstice_id = dst.idx
