import numpy as np


def solve(u, _, conditions):
    x0, t0 = conditions["x_min"], conditions["t_min"]
    xn, tm = conditions["x_max"], conditions["t_max"]
    dt, dx = conditions['dt'], conditions['dx']
    T0 = conditions["T"]

    xs, ts = np.arange(x0, xn, dx), np.arange(t0, tm, dt)
    xl, tl = len(xs), len(ts)

    T = np.zeros([xl, tl])
    for i in range(xl):  # T[x][t]
        T[i][0] = T0(xs[i])

    for j in range(tl - 1):
        for i in range(xl):
            Tj = {
                "i-1": T[__b(i - 1, (0, xl - 1))][j],
                "i": T[i][j],
                "i+1": T[__b(i + 1, (0, xl - 1))][j],
            }

            T[i][j + 1] = __Lh(
                u(ts[j], xs[i], T[i][j]), 0,
                dx, Tj) * 2 * dt + T[i][__b(j - 1, (0, tl - 1))]

    return T


def __b(v, bounds):
    if v < bounds[0]:
        return 0
    if v > bounds[1]:
        return bounds[1]
    return v


def __Lh(u, g, dx, Tj):
    return -u * ((Tj["i+1"] - Tj["i-1"]) / (2 * dx)) + g * ((Tj["i-1"] + Tj["i+1"] - 2 * Tj["i"]) / dx / dx)
