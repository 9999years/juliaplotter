# floor mostly, some assorted stuff
import math
# misc for actually evaluating the julia equations
import cmath
# turning user input into an eval-able formula
import re
# file writing
import codecs
# args? i actually might not need this
import argparse
# warnings
import sys
# making directories, deleting the .ppm after conversion, resolving paths
import os
import random

def warn(*args, **kwargs):
    print('WARNING: ', *args, file=sys.stderr, **kwargs)

# set up arguments
parser = argparse.ArgumentParser(
    description='Render an arbitrary Julia set or a grid of Julia sets'
    'with differing constants c from a user-entered equation.'
)

parser.add_argument('--degree', '-d', position=1, metavar='degree', type=int,
    default='5', help='The polynomial’s maximum degree.')

parser.add_argument('--variables', '-v', position=1, metavar='var', type=str,
    nargs='+', default='x', help='An array of variables to use.')

parser.add_argument('--constants', '-c', position=1, metavar='const', type=str,
    nargs='+', default='x', help='An array of constants to use.')

parser.add_argument('--imaginary', '-i', position=1, action='store_true',
    help='Use imaginary numbers.')

# parse arguments, extract variables
args = parser.parse_args()

degree = args.degree
vars   = args.variables
consts = args.constants
imag   = args.imaginary

top = ""
def term(deg, vars, consts, imag):
    ret = "("
    return ret

def poly(deg, vars, consts, imag):
    ret = ""
    for i in random.randrange(deg):
        ret += term(deg, vars, consts, imag)
