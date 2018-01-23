import numpy as np


def solve(consts, conditions):
    x0, t0 = conditions["z_min"], conditions["t_min"]
    xn, tm = conditions["z_max"], conditions["t_max"]
    dt, dx = conditions['dt'], conditions['dz']

    T0, Tm = conditions["T0"], conditions["Tm"]
    X0, Xn = conditions["X0"], conditions["Xn"]

    D, C, Q = consts["D"], consts["C"], consts["Q"]
    ro, lamda, alpha = consts["ro"], consts["lambda"], consts["alpha"]
    K, E, R = consts["K"], consts["E"], consts["R"]

    xs, ts = np.arange(x0, xn, dx), np.arange(t0, tm, dt)
    xl, tl = len(xs), len(ts)

    X = np.zeros([tl, xl])
    X[0][0] = Xn
    T = np.zeros([tl, xl])
    T[0][0] = Tm
    for i in range(1, xl):
        T[0][i] = T0
        X[0][i] = X0

    for j in range(0, tl - 1):
        # Calc X[n + 1]
        m = np.zeros((xl, xl))
        b = np.zeros(xl)
        for i in range(xl):
            b[i] = X[j][i] / dt
            kappa = D / dx / dx
            w = __w(X[j][i], T[j][i], consts)

            m[i][max(i - 1, 0)] = -kappa
            m[i][min(i + 1, xl - 1)] = -kappa
            m[i][i] = 2 * kappa + 1 / dt + w
        m[0][0] = 1.
        m[0][1] = 0.
        m[-1][-2] = -1.
        m[-1][-1] = 1.
        # print m
        b[0] = 1.
        b[-1] = 0.

        tx = np.linalg.solve(m, b)
        # print np.allclose(np.dot(m, tx), b)
        # print np.dot(m, tx)
        # print tx
        # print np.linalg.solve(m, b)
        # print "--------"
        for i in range(len(tx)):
            X[j + 1][i] = tx[i]

        # Calc T[n + 1]
        m = np.zeros((xl, xl))
        b = np.zeros(xl)

        for i in range(xl):
            # Ld = lamda * dt / ro / C / dx / dx
            w = __w(X[j][i], T[j][i], consts) * X[j][i]
            b[i] = ro * C / dt * T[j][i] - ro * Q * w
            #
            # m[i][max(i - 1, 0)] = -Ld
            # m[i][min(i + 1, xl - 1)] = -Ld
            # m[i][i] = 1 + 2 * Ld
            m[i][max(i - 1, 0)] = -lamda / dx / dx
            m[i][min(i + 1, xl - 1)] = -lamda / dx / dx
            m[i][i] = ro * C / dt + 2 * lamda / dx / dx
        m[0][0] = 1.
        m[0][1] = 0.
        m[-1][-2] = -1.
        m[-1][-1] = 1.

        b[0] = Tm
        b[-1] = 0.

        tx = np.linalg.solve(m, b)
        for i in range(len(tx)):
            T[j + 1][i] = tx[i]

    return X, T


def __w(X, T, consts):
    return consts["K"] * np.power(X, consts["alpha"] - 1.) * np.exp(-consts["E"] / consts["R"] / T)


def __solve_diagonal(m, d):
    n = len(d) - 1
    xs = np.zeros(n + 1)

    als, bls = np.zeros(n + 1), np.zeros(n + 1)
    als[0] = - m[0][1] / m[0][0]
    bls[0] = - d[0] / m[0][0]

    for i in range(1, n):
        a = m[i][i - 1]
        b = m[i][i]
        c = m[i][i + 1]
        als[i] = - c / (a * als[i - 1] + b)
        bls[i] = (d[i] - a * bls[i - 1]) / (a * als[i - 1] + b)

    xs[n] = (d[n] - m[n][n - 1] * m[n - 1][n - 1]) / (m[n][n - 1] * als[n - 1] + m[n][n])

    for i in range(n - 1, -1, -1):
        xs[i] = als[i] * xs[i + 1] + bls[i]
    return xs