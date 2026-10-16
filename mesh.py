import numpy as np


def generate_mesh(nozzle, n):
    x = np.linspace(nozzle.x_start, nozzle.x_end, n)
    dx = nozzle.length / (n - 1)

    #Normalise by the true (continuous) throat area so that A* = 1, even if the throat falls between nodes
    _, area_throat = nozzle.throat()
    area = nozzle.area(x) / area_throat

    return x, dx, area
