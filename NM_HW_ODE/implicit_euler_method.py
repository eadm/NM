import numpy as np


def solve(system, conditions):
    dt = system["dt"]
    t = 0
    x, y, z = conditions
    xs, ys, zs, ts = [x], [y], [z], [t]

    for _ in range(1000):
        f = [
            lambda x1, y1, z1: x + (-system["delta"] * x1 + system["delta"] * y1) * system["dt"],
            lambda x1, y1, z1: y + (-x1 * z1 + system["r"] * x1 - y1) * system["dt"],
            lambda x1, y1, z1: z + (x1 * y1 - system["b"] * z1) * system["dt"]
        ]

        x, y, z = __simple_iterations(f, [x, y, z], 0.01)
        t = t + dt

        xs.append(x)
        ys.append(y)
        zs.append(z)
        ts.append(t)

    return xs, ys, zs, ts


# mb
def __simple_iterations(f, x0, eps):
    x = x0[:]

    dist = 1
    iterations = 10000
    it = 1
    while dist > eps and it < iterations:
        f_v = np.vectorize(lambda fun: fun(x[0], x[1], x[2]))
        x1 = f_v(f)

        dist = np.linalg.norm(x1 - x)
        x = x1
        it += 1

    return x
