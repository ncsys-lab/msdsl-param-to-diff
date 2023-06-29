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
from sympy import Symbol, Derivative
from sympy.integrals import inverse_laplace_transform
#from msdsl import MixedSignalModel, VerilogGenerator, RangeOf
from fixture2differentialeq import *
from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt


THIS_DIR = Path(__file__).resolve().parent
BUILD_DIR = THIS_DIR

class SympyVarGenerator:

    def __init__(self):
        self.vars = {}
        self.derivs = {}
        self.odes = {}

    def symbols(self):
        return self.vars.values()

    def get_odes(self):
        return self.odes.items()

    def add_ode(self,lhs,rhs):
        print("%s' = %s" % (lhs,rhs))
        self.odes[lhs] = rhs

    def derivatives(self):
        for (name,order),expr in self.derivs.items():
            yield order,name,expr

    def highest_order(self,name):
        ords = []
        for n,order in self.derivs.keys():
            if n == name:
                ords.append(order)
        return self.derivs[(name,max(ords))]

    def definevar(self,name):
        assert(not name in self.vars)
        self.vars[name] = Symbol(name=name)
        return self.vars[name]

    def getvar(self,name):
        return self.vars[name]

    def deriv(self,var,order):
        if order == 0:
            return var

        name = var.name
        if not (name,order) in self.derivs:
            for i in range(1,order+1):
                if not (name,i) in self.vars:
                    newvar = self.vars[name] if i - 1 == 0 else self.derivs[(name,i-1)]
                    self.derivs[(name,i)] = Derivative(newvar,t)
        
        return self.derivs[(name,order)]

    def get_assignments(self):
        for (name,i),(var,ddt)  in self.derivs.items():
            yield var == ddt



class GekkoConverter:

    def __init__(self,m):
        self.m = m
        self.vars = {}
        self.conv = {}

    def getvar(self,name):
        return self.vars[name]

    def definevar(self,sym):
        assert(not sym.name in self.vars)
        self.vars[sym.name] = self.m.Var(name=sym.name)
        self.conv[sym] = self.vars[sym.name]
        return self.vars[sym.name]

    def ode(self,lhs,rhs):
        gekko_lhs = self.conv[lhs]
        variables = list(rhs.free_symbols)
        fn = sym.lambdify(variables, rhs)
        args = list(map(lambda v: self.conv[v], variables))
        gekko_rhs = fn(*args)
        return gekko_lhs.dt() == gekko_rhs



def main():
    #test()
    #input("----test finished---")
    #Extract DiffEQ from parameters


    print("constructing diffeq") 
    vargen = SympyVarGenerator()
    polarparam_to_diff_dict(
        circuit_cfg=THIS_DIR / 'circuit.cfg',
        params_yaml=THIS_DIR / 'params.yaml',
        vargen = vargen
    )

    print("---- converting to Gekko ----")
    m = GEKKO()
    gc = GekkoConverter(m)
    for var in vargen.symbols():
        gc.definevar(var)

    for lhs,rhs in vargen.get_odes():
        eqn = gc.ode(lhs,rhs)
        print(eqn)
        m.Equation(eqn)

    #Construct equation. Based on example I found.
    #Can easily be automated, but maybe not best use of time ATM?
    nt = 1000
    TIMESCALE = 1e10
    max_time = 10
    m.time = np.linspace(0,max_time*1.0/TIMESCALE,nt)

    #setting input as step function
    """
    ut = np.zeros(nt)
    ut[2:3] = 1.2
    for i in range(4,5):
        ut[i] = 1.2
    """
    ut = np.sin(2 * pi * TIMESCALE * m.time)
    m.Equation(gc.getvar("du_d1") == m.Param(ut))

    print("-------------")
    input("press any key to solve")
    #m.Equation(ud == u)
    #m.Equation(duddt == ud.dt())
    #m.Equation( Y_dict_t['y'] * y + y.dt() == duddt.dt() + U_dict_t['u\''] * ud.dt() + U_dict_t['u'] / 3e10 * u)

    #Simulation options
    # dynamic sequential
    m.options.IMODE = 7
    m.options.NODES = 3
    m.options.MAX_ITER = 100
    m.solve(disp=False)

    u = gc.getvar("u")
    y = gc.getvar("y")
    #plot the figure!
    plt.figure()
    plt.plot(m.time,u.value,label='u(t)')
    plt.legend()
    plt.xlabel('Time')
    plt.savefig("u.png")
    plt.clf()

    plt.plot(m.time,y.value,label='y(t)')
    plt.legend()
    plt.xlabel('Time')
    plt.savefig("y.png")
    plt.clf()


if __name__ == '__main__':
    main()