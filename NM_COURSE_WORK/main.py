import test_method


print test_method.solve({
    "K": 1.6 * 10e6,
    "E": 8. * 10e4,
    "R": 8.314,
    "alpha": 0.5,
    "Q": 7. * 10e5,
    "ro": 830.,
    "C": 1980.,
    "lambda": 0.13,
    "D": 8. * 10e-12,
    "T0": 293.,
    "Tm": 328.
}, {
    "z_min": 0.,
    "z_max": 10.,
    "t_min": 0.,
    "t_max": 10.,
    "T0": 293.,
    "Tm": 328.,
    "dt": 0.01,
    "dz": 0.01
})
