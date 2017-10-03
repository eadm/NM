def solve(system, conditions):
    dt = system["dt"]

    t = 0
    xi, yi, zi = conditions
    xs, ys, zs, ts = [xi], [yi], [zi], [t]

    for _ in range(1000):
        xi1 = xi + __k(system, "x", xi, yi, zi, dt) * dt
        yi1 = yi + __k(system, "y", xi, yi, zi, dt) * dt
        zi1 = zi + __k(system, "z", xi, yi, zi, dt) * dt

        xi, yi, zi, t = xi1, yi1, zi1, t + dt

        xs.append(xi)
        ys.append(yi)
        zs.append(zi)
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
