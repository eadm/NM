def solve(system, conditions):
    dt = 0.1

    t = 0
    xi, yi, zi = conditions
    xs, ys, zs, ts = [xi], [yi], [zi], [t]

    for _ in range(500):
        xi1 = xi + k(system, "x", xi, yi, zi, dt) * dt
        yi1 = yi + k(system, "y", xi, yi, zi, dt) * dt
        zi1 = zi + k(system, "z", xi, yi, zi, dt) * dt

        xi, yi, zi, t = xi1, yi1, zi1, t + dt

        xs.append(xi)
        ys.append(yi)
        zs.append(zi)
        ts.append(t)

    return xs, ys, zs, ts


def k(system, v, x, y, z, dt):
    return (k0(system, v, x, y, z, dt) + 2 * k1(system, v, x, y, z, dt) +
            2 * k2(system, v, x, y, z, dt) + k3(system, v, x, y, z, dt)) / 6


def k0(system, v, x, y, z, dt):
    return system[v](x, y, z)


def k1(system, v, x, y, z, dt):
    return system[v](x + (dt * k0(system, "x", x, y, z, dt) / 2),
                     y + (dt * k0(system, "y", x, y, z, dt) / 2),
                     z + (dt * k0(system, "z", x, y, z, dt) / 2))


def k2(system, v, x, y, z, dt):
    return system[v](x + (dt * k1(system, "x", x, y, z, dt) / 2),
                     y + (dt * k1(system, "y", x, y, z, dt) / 2),
                     z + (dt * k1(system, "z", x, y, z, dt) / 2))


def k3(system, v, x, y, z, dt):
    return system[v](x + dt * k2(system, "x", x, y, z, dt),
                     y + dt * k2(system, "y", x, y, z, dt),
                     z + dt * k2(system, "z", x, y, z, dt))
