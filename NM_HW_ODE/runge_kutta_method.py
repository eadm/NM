import numpy as np


def solve(system, conditions):
    dt = system["dt"]

    x, y, z = conditions
    xs, ys, zs, ts = [x], [y], [z], [0.]

    for t in np.arange(dt, system["t_max"], dt):
        x1 = x + __k(system, "x", x, y, z, dt) * dt
        y1 = y + __k(system, "y", x, y, z, dt) * dt
        z1 = z + __k(system, "z", x, y, z, dt) * dt

        x, y, z = x1, y1, z1

        xs.append(x)
        ys.append(y)
        zs.append(z)
        ts.append(t)

    return xs, ys, zs, ts


def __k(system, v, x, y, z, dt):
    return (__k0(system, v, x, y, z, dt) + 2 * __k1(system, v, x, y, z, dt) +
            2 * __k2(system, v, x, y, z, dt) + __k3(system, v, x, y, z, dt)) / 6


def __k0(system, v, x, y, z, dt):
    return system[v](x, y, z)


def __k1(system, v, x, y, z, dt):
    return system[v](x + (dt * __k0(system, "x", x, y, z, dt) / 2),
                     y + (dt * __k0(system, "y", x, y, z, dt) / 2),
                     z + (dt * __k0(system, "z", x, y, z, dt) / 2))


def __k2(system, v, x, y, z, dt):
    return system[v](x + (dt * __k1(system, "x", x, y, z, dt) / 2),
                     y + (dt * __k1(system, "y", x, y, z, dt) / 2),
                     z + (dt * __k1(system, "z", x, y, z, dt) / 2))


def __k3(system, v, x, y, z, dt):
    return system[v](x + dt * __k2(system, "x", x, y, z, dt),
                     y + dt * __k2(system, "y", x, y, z, dt),
                     z + dt * __k2(system, "z", x, y, z, dt))
