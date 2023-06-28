#========================================================================
#
#   run.py
#   Janky mish-mash of different examples I found which attempt to convert
#   parameters which convert fixture outputs to a differential equation.
#   
#
#   -w
#
#========================================================================

import yaml
from math import pi
from pathlib import Path

import numpy as np
from sympy.integrals import inverse_laplace_transform
from msdsl import MixedSignalModel, VerilogGenerator, RangeOf
from fixture2differentialeq import *


THIS_DIR = Path(__file__).resolve().parent
BUILD_DIR = THIS_DIR



def main():

    polarparam_to_diff_dict(
        circuit_cfg=THIS_DIR / 'circuit.cfg',
        params_yaml=THIS_DIR / 'params.yaml'
    )




if __name__ == '__main__':
    main()