import numpy as np


def solve(consts, conditions):
    D, C, Q = consts["D"], consts["C"], consts["Q"]
    ro, lamda, alpha = consts["ro"], consts["lambda"], consts["alpha"]
    K, E, R = consts["K"], consts["E"], consts["R"]

    z0, t0 = conditions["z_min"], conditions["t_min"]
    zn, tm = conditions["z_max"], conditions["t_max"]
    dz, dt = conditions['dz'], conditions['dt']
    T0, Tm = conditions["T0"], conditions["Tm"]
    X0, Xn = conditions["X0"], conditions["Xn"]

    zs, ts = np.arange(z0, zn, dz), np.arange(t0, tm, dt)
    zl, tl = len(zs), len(ts)

    T = np.zeros([tl, zl])
    X = np.zeros([tl, zl])

    def W(x, t):
        return -K * (x ** alpha) * np.math.exp(-E / R / t)

    for i in range(tl):  # T[t][x]
        X[i][0] = X0
        X[i][zl - 1] = Xn
        T[i][0] = Tm
        T[i][zl - 1] = T0

    for j in range(1, zl):
        X[0][j] = 1
        T[0][j] = T0

    for j in range(tl - 1):
        for i in range(1, zl - 1):
            Xj = {
                "i-1": X[j][__b(i - 1, (0, zl - 1))],
                "i": X[j][i],
                "i+1": X[j][__b(i + 1, (0, zl - 1))],
            }

            X[j + 1][i] = (W(X[j][i], T[j][i]) + D * _double_d(Xj, dz)) * dt + X[j][i]

            Tj = {
                "i-1": T[j][__b(i - 1, (0, zl - 1))],
                "i": T[j][i],
                "i+1": T[j][__b(i + 1, (0, zl - 1))],
            }

            T[j + 1][i] = ((-Q / C) * W(X[j][i], T[j][i]) + (lamda / ro / C) * _double_d(Tj, dz)) * dt + T[j][i]

    return X, T


def __b(v, bounds):
    if v < bounds[0]:
        return 0
    if v > bounds[1]:
        return bounds[1]
    return v


def _double_d(Uj, dz):
    return (Uj["i-1"] + Uj["i+1"] - 2 * Uj["i"]) / dz / dz
