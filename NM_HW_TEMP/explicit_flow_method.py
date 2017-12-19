import numpy as np


def solve(u, g, conditions):
    x0, t0 = conditions["x_min"], conditions["t_min"]
    xn, tm = conditions["x_max"], conditions["t_max"]
    dt, dx = conditions['dt'], conditions['dx']
    T0 = conditions["T"]

    xs, ts = [], []
    for i in np.arange(x0, xn, dx):
        xs.append(i)
    for j in np.arange(t0, tm, dt):
        ts.append(j)
    xs, ts = np.array(xs), np.array(ts)
    xl, tl = len(xs), len(ts)

    T = np.zeros([xl, tl])
    for i in range(xl):  # T[x][t]
        T[i][0] = T0(xs[i])

    for j in range(0, tl - 1):
        for i in range(xl):
            Tj = {
                "i-1": T[__b(i - 1, (0, xl - 1))][j],
                "i": T[i][j],
                "i+1": T[__b(i + 1, (0, xl - 1))][j],
            }

            T[i][j + 1] = __Lh(
                u(ts[j], xs[i], T[i][j]),
                g(ts[j], xs[i], T[i][j]),
                dx, Tj) * dt + T[i][j]

    return T


def __b(v, bounds):
    if v < bounds[0]:
        return 0
    if v > bounds[1]:
        return bounds[1]
    return v


def __u(u, t, x, T):
    return u(t, x, T)


def __g(g, t, x, T):
    return g(t, x, T)


def __Lh(u, g, dx, Tj):
    return -u * ((Tj["i"] - Tj["i-1"]) / dx) + g * ((Tj["i-1"] + Tj["i+1"] - 2 * Tj["i"]) / dx / dx)
