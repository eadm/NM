# coding=utf-8
from mpl_toolkits.mplot3d import Axes3D
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


py.plot(ts, xs, label="x(t)")
py.plot(ts, ys, label="y(t)")
py.plot(ts, zs, label="z(t)")

fig = py.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(xs, ys, zs, label="x,y,z")

py.show()
