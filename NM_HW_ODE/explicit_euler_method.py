# import math
# import numpy as __np


def solve(system, conditions):
    dt = 0.1

    t = 0

    x, y, z = conditions
    xs, ys, zs, ts = [x], [y], [z], [t]

    for _ in range(500):
        x1 = x + system["x"](x, y, z) * dt
        y1 = y + system["y"](x, y, z) * dt
        z1 = z + system["z"](x, y, z) * dt

        x, y, z, t = x1, y1, z1, t + dt

        xs.append(x)
        ys.append(y)
        zs.append(z)
        ts.append(t)

    return xs, ys, zs, ts
