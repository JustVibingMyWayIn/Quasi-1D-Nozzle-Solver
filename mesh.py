#mesh generation mesh.py
import numpy as np

def generate_mesh(n, length):
    x = np.linspace(0.0, length, n)
    dx = length / (n - 1)

    area = 1 + 2.2 * (x - 1.5) ** 2
    area /= np.min(area)

    return x, dx, area