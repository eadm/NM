import numpy as np


def solve(u, g, conditions):
    x0, t0 = conditions["x_min"], conditions["t_min"]
    xn, tm = conditions["x_max"], conditions["t_max"]
    dt, dx = conditions['dt'], conditions['dx']
    T0 = conditions["T"]

    xs, ts = np.arange(x0, xn, dx), np.arange(t0, tm, dt)
    xl, tl = len(xs), len(ts)

    T = np.zeros([xl, tl])
    for i in range(xl):  # T[x][t]
        T[i][0] = T0(xs[i])

    for j in range(0, tl - 1):
        m = np.zeros((xl, xl))
        b = np.zeros(xl)
        for i in range(xl):
            s = u(ts[j], xs[i], T[i][j]) * dt / dx
            r = g(ts[j], xs[i], T[i][j]) * dt / dx / dx

            b[i] = T[i][j]
            m[i][max(i - 1, 0)] = -(s + r)
            m[i][min(i + 1, xl - 1)] = -r
            m[i][i] = (1 + s + 2 * r)

        ts = __solve_diagonal(m, b)
        for i in range(len(ts)):
            T[i][j + 1] = ts[i]

    return T


def __solve_diagonal(m, b):

    return []
