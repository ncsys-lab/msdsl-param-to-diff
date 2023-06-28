#========================================================================
#
#   fixture2differentialeq.py
#   Functions which take parameters from fixture and output a differential
#   equation which describes the sytstem
#
#   -w
#
#========================================================================

import yaml
from math import pi
from pathlib import Path

import numpy as np
from sympy.integrals import inverse_laplace_transform
import sympy as sym
from sympy.abc import s,t,x,y,z

#Laplace and time domain dictionary... add more as needed - can probably chatGPT this
laplace_time_dict = {
    1 : "x",
    s : "x'",
    s**2 : "x''"
}

def read_yaml( f ):
    with open(f, 'r') as stream:
        return yaml.safe_load(stream)



def polarparam_to_coeff_dict( circuit_cfg, params_yaml ):
        

        # read the two yaml files
        circuit = read_yaml( circuit_cfg )
        params  = read_yaml( params_yaml )

        #This lambda indexes into the different parameters of the analog block
        read_param = lambda name: params['test1'][name][0]['coef']['1']

        dcgain = read_param( 'dcgain' )

        #Take the poles from the file
        #the msdsl demo uses some sort of timestep-based solver, here I just
        #solve the entire thing using sympy
                
        p1 = read_param( 'fp1' )
        p2 = read_param( 'fp2' )
        z1 = read_param( 'fz1' )

        #Initialize the transfer function, need to manually enter this for the time being.
        G = ( s - z1 ) / ( ( s - p1 ) * ( s - p2 ) )

        #Simplify G, just in case
        G = sym.simplify(G)

        # G = Y(s)/U(s). Isolate both sides
        Y_side, U_side = sym.fraction(G)

        #Expand Factors
        Y_side = sym.expand_mul(Y_side)
        U_side = sym.expand_mul(U_side)
        print(Y_side)
        print(U_side)

        #create dictionary of coefficients
        Y_dict = Y_side.as_coefficients_dict(s)
        U_dict = U_side.as_coefficients_dict(s)
        print(Y_dict)
        print(U_dict)

        #return for further processing by other means.
        return (Y_dict, U_dict)

#equation_tuple is ( Y_dict, U_dict )
def laplace_domain_dict_to_diff_dict( equation_tuple ):
    Y_dict = equation_tuple[0]
    U_dict = equation_tuple[1]

    #python doesn't like when you try to change dictionary mid-loop
    Y_dict_t = {}
    U_dict_t = {}

    #output side
    for term in Y_dict:
        #replace each key with time-domain equivalent
        Y_dict_t[laplace_time_dict[term].replace("x","y")] = Y_dict[term]
    
    #input side
    for term in U_dict:

        U_dict_t[laplace_time_dict[term].replace("x","u")] = U_dict[term]


    print(( Y_dict_t, U_dict_t ))
    return ( Y_dict_t, U_dict_t )


def polarparam_to_diff_dict( circuit_cfg, params_yaml ):
    coeff_tuple = polarparam_to_coeff_dict( circuit_cfg, params_yaml)
    return laplace_domain_dict_to_diff_dict(coeff_tuple)
