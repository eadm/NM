import numpy as np
import explicit_flow_method
import pylab as pl
import explicit_counter_flow_method

cs = explicit_flow_method.solve(
    lambda t, x, T: 1.3,
    lambda t, x, T: 0.,
    {
        "t_min": 0.,
        "t_max": 20.,
        "dt": 0.2,
        "dx": 1.,
        "T": lambda x: x ** 2,
        "x_min": 0.,
        "x_max": 20.
    }
)

m_x, m_y = np.mgrid[slice(0., 20., 1.), slice(0., 20., 0.2)]


pl.pcolormesh(m_x, m_y, cs)
pl.show()
print cs
