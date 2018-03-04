import numpy as np
import pylab as pl
import subprocess

import c_wrapper
import parser

consts = {
    "D": 8e-12,
    "C": 1980.,
    "Q": 7e5,

    "ro": 830.,
    "lambda": 0.13,
    "alpha": 1.,

    "K": 1.6e6,
    "E": 8e4,
    "R": 8.314,
}

conditions = {
    "z_min": 0.,
    "z_max": 0.03,
    "dz": 0.0001,

    "t_min": 0.,
    "t_max": 500.,
    "dt": 0.01,

    "X0": 1.,
    "Xn": 0.,

    "T0": 293.,
    "Tm": 293. + (consts["Q"] / consts["C"])
}

print "Tm = " + str(conditions["Tm"]) + " K"

kappa = consts["lambda"] / consts["ro"] / consts["C"]


# consts["D"] = kappa


def v():
    f = 2 * consts["K"] * consts["lambda"] / consts["Q"] / consts["ro"] / (conditions["Tm"] - conditions["T0"])
    s = (consts["R"] * (conditions["Tm"] ** 2) / consts["E"]) ** 2
    t = np.exp(-consts["E"] / consts["R"] / conditions["Tm"])
    return np.sqrt(f * s * t)


print "U = " + str(v()) + " m/sec"
print "U = " + str(v() * 6000) + " mm/min"

print "kappa = " + str(kappa)
sigma_h = kappa / v()
print "sigma_h = " + str(sigma_h)
sigma_D = consts["D"] / v()
print "sigma_D = " + str(sigma_D)
betta = consts["R"] * conditions["Tm"] / consts["E"]
print "betta = " + str(betta)
sigma_r = sigma_h * betta
print "sigma_r = " + str(sigma_r)

l_args = c_wrapper.get_args('./c_lib/nm_course_work', consts, conditions)

p = subprocess.Popen(l_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output = p.communicate()[0]
print output

X = parser.parse_matrix('X.out')
T = parser.parse_matrix('T.out')

print X
print T

m_z, m_y = np.mgrid[
    slice(conditions["t_min"], conditions["t_max"], conditions["dt"]),
    slice(conditions["z_min"], conditions["z_max"], conditions["dz"])]

pl.pcolormesh(m_z, m_y, T)
pl.show()
