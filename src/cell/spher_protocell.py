import numpy as np
from src.cell.protocell import Protocell

#
#   Spherical protocell class
#

class SphProtocell(Protocell):
    def __init__(self, idx, radius, n_shells, contact_to_mineral_surface=False):
        super().__init__(idx=idx, contact_to_mineral_surface=contact_to_mineral_surface)
        # structural data
        self.geometry = "spherical"
        self.radius = radius
        self.n_shells = n_shells
        # set volume
        self.volume = 4./3 * np.pi * self.radius ** 3
        # set shell grid
        self.sh_edges = None
        self.sh_centers = None
        self.shell_volumes = None
        self.build_radial_grid()
    def build_radial_grid(self):
        """
        Partition protocell volume into n_cells spherical shells.
        """
        # shell edges and centers
        self.sh_edges = np.linspace(0.0, self.radius, self.n_shells + 1)
        self.sh_centers = 0.5 * (self.sh_edges[:-1] + self.sh_edges[1:])
        # shell volumes
        self.shell_volumes = (4.0 / 3.0) * np.pi * (
            self.sh_edges[1:]**3 - self.sh_edges[:-1]**3
        )
    def sample_relative_metabolite_positions(count, rng):
        """
        Sample particle positions uniformly inside a sphere centered at the origin.
        Returns positions relative to protocell center.
        """
        if rng is None:
            rng = np.random.default_rng()
        print(rng)
        # protocell radius
        self.radius = (3.0 * self.volume / (4.0 * np.pi)) ** (1./3.)
        # define position inside sphere
        u = rng.random(n_particles)
        cos_theta = rng.uniform(-1.0, 1.0, n_particles)
        phi = rng.uniform(0.0, 2.0 * np.pi, n_particles)
        r = radius * u ** (1./3)
        sin_theta = np.sqrt(1.0 - cos_theta**2)
        # cartesian
        x = r * sin_theta * np.cos(phi)
        y = r * sin_theta * np.sin(phi)
        z = r * cos_theta
        # molecule distr.
        rho = np.column_stack((x, y, z))
        return rho