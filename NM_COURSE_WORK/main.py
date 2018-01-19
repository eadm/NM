import test_method
import numpy as np
import pylab as pl

consts = {
    "K": 1.6 * 10e6,
    "E": 8. * 10e4,
    "R": 8.314,
    "alpha": 0.5,
    "Q": 7. * 10e5,
    "ro": 830.,
    "C": 1980.,
    "lambda": 0.13,
    "D": 8. * 10e-12
}

conditions = {
    "z_min": 0.,
    "z_max": 10.,
    "t_min": 0.,
    "t_max": 10.,
    "T0": 293.,
    "Tm": 328.,
    "X0": 0.,
    "Xn": 1.,
    "dt": 0.01,
    "dz": 0.01
}

X, T = test_method.solve(consts, conditions)

print X, T

m_z, m_y = np.mgrid[
    slice(conditions["z_min"], conditions["z_max"], conditions["dz"]),
    slice(conditions["t_min"], conditions["t_max"], conditions["dt"])]

# print T
pl.pcolormesh(m_z, m_y, X)
pl.show()