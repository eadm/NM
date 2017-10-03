import numpy as np


def solve(system, conditions, dt):
    t = 0
    x, y, z = conditions
    xs, ys, zs, ts = [x], [y], [z], [t]

    for _ in range(1000):
        F = [
            lambda x1, y1, z1: x - (system["delta"] * y1 + x) / (1 + system["delta"] * system["dt"]),
            lambda x1, y1, z1: y - (y + system["r"] * x1 * system["dt"] - x1 * z1 * system["dt"]) / (1 + system["dt"]),
            lambda x1, y1, z1: z - (z + x1 * y1 * system["b"] * system["dt"]) / (1 + system["b"] * system["dt"])
        ]

        FM = [
            [
                lambda x1, y1, z1: 1,
                lambda x1, y1, z1: -(system["delta"] * system["dt"]) / (1 + system["delta"] * system["dt"]),
                lambda x1, y1, z1: 0],
            [
                lambda x1, y1, z1: -(system["r"] - z1) * system["dt"]/(1 + system["dt"]),
                lambda x1, y1, z1: 1,
                lambda x1, y1, z1: x1 * system["dt"] / (1 + system["dt"])
            ],
            [
                lambda x1, y1, z1: -y1 * system["b"] * system["dt"] / (1 + system["b"] * system["dt"]),
                lambda x1, y1, z1: -x1 * system["b"] * system["dt"] / (1 + system["b"] * system["dt"]),
                lambda x1, y1, z1: 1
            ]
        ]

        x, y, z = __newton_method(F, FM, [x, y, z], 0.1)
        print [x, y, z]
        t = t + dt

        xs.append(x)
        ys.append(y)
        zs.append(z)
        ts.append(t)


    return xs, ys, zs, ts


def __newton_method(F, FM, x0, eps):
    x = x0[:]

    # x = x -
    print "x:"
    print x

    dist = 1
    while dist > eps:
        f_v = np.vectorize(lambda f: f(x[0], x[1], x[2]))
        W = f_v(FM)
        F1 = f_v(F)

        x1 = x - np.linalg.inv(W).dot(F1)

        dist = np.linalg.norm(x1 - x)
        x = x1
        print x

    return x
