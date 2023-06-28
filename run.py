#========================================================================
#
#   run.py
#   Extracts parameters and simulates differential equation
#   uses GEKKO
#
#========================================================================

import yaml
from math import pi
from pathlib import Path

import numpy as np
from sympy.integrals import inverse_laplace_transform
from msdsl import MixedSignalModel, VerilogGenerator, RangeOf
from fixture2differentialeq import *
from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt


THIS_DIR = Path(__file__).resolve().parent
BUILD_DIR = THIS_DIR

m = GEKKO()
nt = 1000
TIMESCALE = 1e10
m.time = np.linspace(0,1/TIMESCALE,nt)

#setting input as step function
"""
ut = np.zeros(nt)
ut[2:3] = 1.2
for i in range(4,5):
    ut[i] = 1.2
"""
ut = np.sin(2 * pi * TIMESCALE * m.time)

def main():
    #Extract DiffEQ from parameters
    ( Y_dict_t, U_dict_t ) = polarparam_to_diff_dict(
        circuit_cfg=THIS_DIR / 'circuit.cfg',
        params_yaml=THIS_DIR / 'params.yaml'
    )
    
    #Construct equation. Based on example I found.
    #Can easily be automated, but maybe not best use of time ATM?
    u = m.Param(value = ut)
    ud = m.Var()
    duddt = m.Var()

    y = m.Var()
    
    m.Equation(ud == u)
    m.Equation(duddt == ud.dt())
    m.Equation( Y_dict_t['y'] * y + y.dt() == duddt.dt() + U_dict_t['u\''] * ud.dt() + U_dict_t['u'] / 3e10 * u)

    #Simulation options
    m.options.IMODE = 7
    m.options.NODES = 3
    m.options.MAX_ITER = 100
    m.solve(disp=True)

    #plot the figure!
    plt.figure()
    plt.plot(m.time,u.value,label='u(t)')
    plt.plot(m.time,y.value,label='y(t)')
    plt.legend()
    plt.xlabel('Time')
    plt.show()


if __name__ == '__main__':
    main()