# coding=utf-8
import pylab as py
import explicit_euler_method

delta = 10.
b = 8. / 3
r = 4

system = {
    "x": lambda x, y, z: -delta * x + delta * y,
    "y": lambda x, y, z: -x * z + r * x - y,
    "z": lambda x, y, z: x * y - b * z
}

# начальные условия (x, y, z)
conditions = (1, 2, 3)

xs, ys, zs, ts = explicit_euler_method.solve(system, conditions)


py.plot(xs, ts)
py.show()
