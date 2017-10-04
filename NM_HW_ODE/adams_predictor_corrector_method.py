import numpy as np
import runge_kutta_method


def solve(system, conditions):
    dt = system["dt"]
    t_max = system["t_max"]

    system["t_max"] = 3 * dt
    xs, ys, zs, ts = runge_kutta_method.solve(system, conditions)

    system["t_max"] = t_max
    t = ts[-1] + dt

    for t in np.arange(t, system["t_max"], dt):
        xp, yp, zp = __predictor(system, xs, ys, zs, dt)
        x, y, z = __corrector(system, xs, ys, zs, xp, yp, zp, dt)

        xs.append(x)
        ys.append(y)
        zs.append(z)
        ts.append(t)

    return xs, ys, zs, ts


def __predictor(system, xs, ys, zs, dt):
    xp = xs[-1] + (23 * system["x"](xs[-1], ys[-1], zs[-1]) - 16 * system["x"](xs[-2], ys[-2], zs[-2])
                   + 5 * system["x"](xs[-3], ys[-3], zs[-3])) * dt / 12
    yp = ys[-1] + (23 * system["y"](xs[-1], ys[-1], zs[-1]) - 16 * system["y"](xs[-2], ys[-2], zs[-2])
                   + 5 * system["y"](xs[-3], ys[-3], zs[-3])) * dt / 12
    zp = ys[-1] + (23 * system["z"](xs[-1], ys[-1], zs[-1]) - 16 * system["z"](xs[-2], ys[-2], zs[-2])
                   + 5 * system["z"](xs[-3], ys[-3], zs[-3])) * dt / 12
    return xp, yp, zp


def __corrector(system, xs, ys, zs, xp, yp, zp, dt):
    x1 = xs[-1] + (5 * system["x"](xp, yp, zp) + 8 * system["x"](xs[-1], ys[-1], zs[-1])
                   - system["x"](xs[-2], ys[-2], zs[-2])) * dt / 12
    y1 = ys[-1] + (5 * system["y"](xp, yp, zp) + 8 * system["y"](xs[-1], ys[-1], zs[-1])
                   - system["y"](xs[-2], ys[-2], zs[-2])) * dt / 12
    z1 = zs[-1] + (5 * system["z"](xp, yp, zp) + 8 * system["z"](xs[-1], ys[-1], zs[-1])
                   - system["z"](xs[-2], ys[-2], zs[-2])) * dt / 12

    return x1, y1, z1
