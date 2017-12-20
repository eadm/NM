import numpy as np


def get_flash_light():
    u = lambda t, x, T: 0.
    g = lambda t, x, T: 0.8
    conditions = {
        "t_min": 0.,
        "t_max": 50.,
        "dt": 0.05,
        "dx": 0.5,
        "T": lambda x: np.math.exp(-abs(x / 10)) * 10,
        "x_min": -100.,
        "x_max": 100.
    }
    return u, g, conditions


def get_steps():
    u = lambda t, x, T: 0
    g = lambda t, x, T: 0.6
    conditions = {
        "t_min": 0.,
        "t_max": 200.,
        "dt": 0.2,
        "dx": 0.5,
        "T": lambda x: x > 0,
        "x_min": -200.,
        "x_max": 200.
    }
    return u, g, conditions


def get_waves():
    u = lambda t, x, T: 5.
    g = lambda t, x, T: 0.
    conditions = {
        "t_min": 0.,
        "t_max": 50.,
        "dt": 0.05,
        "dx": 0.5,
        "T": lambda x: (np.math.cos(x / 20) + 1) * 50,
        "x_min": -100.,
        "x_max": 100.
    }
    return u, g, conditions


def get_rhombus():
    u = lambda t, x, T: 10.
    g = lambda t, x, T: 0.
    conditions = {
        "t_min": 0.,
        "t_max": 50.,
        "dt": 0.05,
        "dx": 0.5,
        "T": lambda x: x > 0,
        "x_min": -100.,
        "x_max": 100.
    }
    return u, g, conditions
