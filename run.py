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
from scipy.integrate import odeint


THIS_DIR = Path(__file__).resolve().parent
BUILD_DIR = THIS_DIR

class SympyVarGenerator:

    def __init__(self):
        self.vars = {}
        self.vars_init_value = {}
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

    def definevar(self,name,init_value=0):
        assert(not name in self.vars)
        self.vars[name] = Symbol(name=name)
        if(init_value != 0):
            self.vars_init_value[name] = init_value
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

    def definevar(self,sym,init_value = 0):
        assert(not sym.name in self.vars)
        self.vars[sym.name] = self.m.Var(name=sym.name, value = init_value)
        self.conv[sym] = self.vars[sym.name]
        return self.vars[sym.name]

    def ode(self,lhs,rhs):
        gekko_lhs = self.conv[lhs]
        variables = list(rhs.free_symbols)
        fn = sym.lambdify(variables, rhs)
        args = list(map(lambda v: self.conv[v], variables))
        gekko_rhs = fn(*args)
        return gekko_lhs.dt() == gekko_rhs

class NumpyODEConverter:
    def __init__(self):
        self.diffeqs = {}
        self.varorder = {}
        self.sim_varorder = []

    def ode(self,lhs,rhs):
        varname = lhs.name
        variables = list(rhs.free_symbols)
        self.diffeqs[varname] = sym.lambdify(variables, rhs)
        self.varorder[varname] = list(map(lambda v: v.name, variables))
        self.sim_varorder.append(varname)

    def compute(self, value_dict):
        result_dict = {}
        for var, diff in self.diffeqs.items():
            variable_order = self.varorder[var]
            variable_list = list(map(lambda v: value_dict[v], variable_order)) 
            result_dict[var] = self.diffeqs[var](*variable_list)
        return result_dict
    
    def scipy_ddt(self, y, t):
        y_dict = dict(map(lambda v: ( self.sim_varorder[v[0]], v[1] ), enumerate(y)))
        y_dict['du_d5'] = 0
        
        if(t > 7e-6 / 2):
            y_dict['u'] = 0
        result_dict = self.compute(y_dict)
        if(t > 7e-6 / 2):
            result_dict['u'] = 0
        return list(map(lambda v: result_dict[v], self.sim_varorder))
    
    def scipy_init(self, i):
        return list(map(lambda v: i[v], self.sim_varorder))
    
    



def main():
    #test()
    #input("----test finished---")
    #Extract DiffEQ from parameters


    print("constructing diffeq") 
    vargen = SympyVarGenerator()
    polarparam_to_diff_dict(
        circuit_cfg=THIS_DIR / 'circuit.cfg',
        params_yaml=THIS_DIR / 'regression_results_manypole.yaml',
        vargen = vargen,
        default_u_value = 1,
        default_y_value = 3.3
    )


    ode_converter = NumpyODEConverter()

    for lhs,rhs in vargen.get_odes():
        ode_converter.ode(lhs,rhs)
        print(lhs)
        print(rhs)

    init_cond = {"u" : 1.6,
                 "y" : 3.3,
                 "du_d1": 0,
                 "du_d2": 0,
                 "du_d3": 0,
                 "du_d4": 0,
                 "dy_d1": 0,
                 "dy_d2": 0,
                 "dy_d3": 0,
                 "dy_d4": 0,
                 "dy_d5": 0,
                 "dy_d6": 0
}

    tmax = 7e-6

    t = np.linspace(0,tmax,1000)

    def ddt(y,t,diffslv):
        result = diffslv.scipy_ddt(y,t)

        return result
    
    sol = odeint(ddt, ode_converter.scipy_init(init_cond), t, args=(ode_converter,))

    for i, var in enumerate(ode_converter.sim_varorder):
        if('dy' in var or 'du' in var):
            continue
        if(var == 'u'):
            print(sol[:, i])
        plt.plot(t, sol[:, i], label=var)
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()

    raise Exception("Stop")

    u_bias = -1
    print("---- converting to Gekko ----")
    m = GEKKO()
    gc = GekkoConverter(m)
    for var in vargen.symbols():
        print(var)
        if(var == Symbol(name='u')):
            gc.definevar(var,init_value=u_bias)
            print('in u loop')
        elif(var == Symbol(name='y')):
            gc.definevar(var,init_value=3.3)
            print('in u loop')
        elif(var == Symbol(name='dy_d1')):
            gc.definevar(var,init_value=1)
            print('in u loop')
        else:
            gc.definevar(var)

    for lhs,rhs in vargen.get_odes():
        eqn = gc.ode(lhs,rhs)
        print(eqn)
        m.Equation(eqn)

    #Construct equation. Based on example I found.
    #Can easily be automated, but maybe not best use of time ATM?
    nt = 10000
    TIMESCALE = 1e6
    max_time = 10
    m.time = np.linspace(0,max_time*1.0/TIMESCALE,nt)

    #setting input as step function
    """
    ut = np.zeros(nt)
    ut[2:3] = 1.2
    for i in range(4,5):
        ut[i] = 1.2
    """
    ut = u_bias
    print(gc.getvar("u"))
    m.Equation(gc.getvar("u") == m.Param(ut))

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
    plt.show()
    plt.savefig("u.png")
    plt.clf()
    
    plt.plot(m.time,y.value,label='y(t)')
    plt.legend()
    plt.xlabel('Time')
    plt.savefig("y.png")
    plt.show()
    plt.clf()



if __name__ == '__main__':
    main()