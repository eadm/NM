import test_method
import numpy as np
import pylab as pl

import implicit

consts = {
    "K": 1.6e6,
    "E": 8.e4,
    "R": 8.314,
    "alpha": 2.,
    "Q": 7.e5,
    "ro": 830.,
    "C": 1980.,
    "lambda": 0.13,
    "D": 8.e-12
}

conditions = {
    "z_min": 0.,
    "z_max": 0.1,
    "t_min": 0.,
    "t_max": 300.,
    "T0": 293.,
    "Tm": 293. + (consts["Q"] / consts["C"]),
    "X0": 0.,
    "Xn": 1.,
    "dt": 0.3,
    "dz": 0.00001
}

print "Tm = " + str(conditions["Tm"]) + " K"


def v():
    f = 2 * consts["K"] * consts["lambda"] / consts["Q"] / consts["ro"] / (conditions["Tm"] - conditions["T0"])
    s = (consts["R"] * (conditions["Tm"] ** 2) / consts["E"]) ** 2
    t = np.exp(-consts["E"] / consts["R"] / conditions["Tm"])
    return np.sqrt(f * s * t)


print "U = " + str(v()) + " m/sec"
print "U = " + str(v() * 6000) + " mm/min"

kappa = consts["lambda"] / consts["ro"] / consts["C"]

print "kappa = " + str(kappa)

sigma_h = kappa / v()
print "sigma_h = " + str(sigma_h)
sigma_D = consts["D"] / v()
print "sigma_D = " + str(sigma_D)
betta = consts["R"] * conditions["Tm"] / consts["E"]
print "betta = " + str(betta)
sigma_r = sigma_h * betta
print "sigma_r = " + str(sigma_r)

X, T = implicit.solve(consts, conditions)

print X
print T

m_z, m_y = np.mgrid[
    slice(conditions["t_min"], conditions["t_max"], conditions["dt"]),
    slice(conditions["z_min"], conditions["z_max"], conditions["dz"])]

# print T
pl.pcolormesh(m_z, m_y, X)
pl.show()
