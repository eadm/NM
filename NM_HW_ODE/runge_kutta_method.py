def solve(system, conditions):
    dt = 0.1

    t = 0
    xi, yi, zi = conditions
    xs, ys, zs, ts = [xi], [yi], [zi], [t]

    k0 = {
        "x": lambda x, y, z: system["x"](x, y, z),
        "y": lambda x, y, z: system["y"](x, y, z),
        "z": lambda x, y, z: system["z"](x, y, z)
    }
    k1 = {
        "x": lambda x, y, z: system["x"](x + (dt * k0["x"](x, y, z) / 2),
                                         y + (dt * k0["y"](x, y, z) / 2),
                                         z + (dt * k0["z"](x, y, z) / 2)
                                         ),
        "y": lambda x, y, z: system["y"](x + (dt * k0["x"](x, y, z) / 2),
                                         y + (dt * k0["y"](x, y, z) / 2),
                                         z + (dt * k0["z"](x, y, z) / 2)
                                         ),
        "z": lambda x, y, z: system["z"](x + (dt * k0["x"](x, y, z) / 2),
                                         y + (dt * k0["y"](x, y, z) / 2),
                                         z + (dt * k0["z"](x, y, z) / 2)
                                         )
    }
    k2 = {
        "x": lambda x, y, z: system["x"](x + (dt * k1["x"](x, y, z) / 2),
                                         y + (dt * k1["y"](x, y, z) / 2),
                                         z + (dt * k1["z"](x, y, z) / 2)
                                         ),
        "y": lambda x, y, z: system["y"](x + (dt * k1["x"](x, y, z) / 2),
                                         y + (dt * k1["y"](x, y, z) / 2),
                                         z + (dt * k1["z"](x, y, z) / 2)
                                         ),
        "z": lambda x, y, z: system["z"](x + (dt * k1["x"](x, y, z) / 2),
                                         y + (dt * k1["y"](x, y, z) / 2),
                                         z + (dt * k1["z"](x, y, z) / 2)
                                         )
    }
    k3 = {
        "x": lambda x, y, z: system["x"](x + dt * k1["x"](x, y, z),
                                         y + dt * k1["y"](x, y, z),
                                         z + dt * k1["z"](x, y, z)
                                         ),
        "y": lambda x, y, z: system["y"](x + dt * k1["x"](x, y, z),
                                         y + dt * k1["y"](x, y, z),
                                         z + dt * k1["z"](x, y, z)
                                         ),
        "z": lambda x, y, z: system["z"](x + dt * k1["x"](x, y, z),
                                         y + dt * k1["y"](x, y, z),
                                         z + dt * k1["z"](x, y, z)
                                         )
    }

    k = lambda v, x, y, z: (k0[v](x, y, z) + 2 * k1[v](x, y, z) + 2 * k2[v](xi, yi, zi) + k3[v](x, y, z)) / 6

    for _ in range(500):
        xi1 = xi + k("x", xi, yi, zi) * dt
        yi1 = yi + k("y", xi, yi, zi) * dt
        zi1 = zi + k("z", xi, yi, zi) * dt

        xi, yi, zi, t = xi1, yi1, zi1, t + dt

        xs.append(xi)
        ys.append(yi)
        zs.append(zi)
        ts.append(t)

    return xs, ys, zs, ts
