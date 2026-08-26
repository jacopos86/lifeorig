from dataclasses import dataclass
import numpy as np


@dataclass
class SpatialGrid1D:
    n_cells: int
    cell_edges: object
    cell_centers: object
    cell_widths: object

#
#   define uniform grid
#

def uniform_spatial_grid_1d(length, n_cells):
    if n_cells <= 0:
        raise ValueError("n_cells must be a positive integer")
    # get grid
    length_units = length.units
    length_values = np.linspace(0.0, length.magnitude, n_cells + 1)
    cell_edges = length_values * length_units
    cell_centers = 0.5 * (cell_edges[:-1] + cell_edges[1:])
    cell_widths = cell_edges[1:] - cell_edges[:-1]
    # return spatial grid
    return SpatialGrid1D(
        n_cells=n_cells,
        cell_edges=cell_edges,
        cell_centers=cell_centers,
        cell_widths=cell_widths,
    )
