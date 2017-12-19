import numpy as np


def get_flash_light():
    u = lambda t, x, T: 0.
    g = lambda t, x, T: 0.6
    conditions = {
        "t_min": 0.,
        "t_max": 20.,
        "dt": 0.2,
        "dx": 0.5,
        "T": lambda x: np.math.exp(-abs(x)),
        "x_min": -20.,
        "x_max": 20.
    }
    return u, g, conditions
