import numpy as np


def solve(system, conditions):
    dt = system["dt"]

    x, y, z = conditions
    xs, ys, zs, ts = [x], [y], [z], [0.]

    for t in np.arange(dt, system["t_max"], dt):
        x1 = x + system["x"](x, y, z) * dt
        y1 = y + system["y"](x, y, z) * dt
        z1 = z + system["z"](x, y, z) * dt

        x, y, z = x1, y1, z1

        xs.append(x)
        ys.append(y)
        zs.append(z)
        ts.append(t)

    return xs, ys, zs, ts
