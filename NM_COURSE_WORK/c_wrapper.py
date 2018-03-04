def get_args(c_path, consts, conditions):
    return [
        c_path,
        str(conditions['z_min']),
        str(conditions['z_max']),
        str(conditions['dz']),
        str(conditions['t_min']),
        str(conditions['t_max']),
        str(conditions['dt']),
        str(conditions['X0']),
        str(conditions['Xn']),
        str(conditions['T0']),
        str(conditions['Tm']),
        str(consts['D']),
        str(consts['C']),
        str(consts['Q']),
        str(consts['ro']),
        str(consts['lambda']),
        str(consts['alpha']),
        str(consts['K']),
        str(consts['E']),
        str(consts['R'])
    ]
