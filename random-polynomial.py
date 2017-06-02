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
# random rationals
import random

def warn(*args, **kwargs):
    print('WARNING: ', *args, file=sys.stderr, **kwargs)

def sign(x):
    return math.copysign(1, x)

math.sign = sign

def signstr(num):
    if math.sign(num) == -1:
        return '-'
    else:
        return '+'

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

def get_coef(imag):
    max_coef = 20
    coef = (rand.random() - 0.5) * max_coef
    if imag:
        coef += get_coef(False) * 1j
    return coef

def term(deg, vars, consts, imag):
    ret = "("
    deg = rand.randint(0, deg)
    var = rand.choice(vars)
    coef = get_coef(imag)
    ret += (f'({rand.choice(vars)}^{deg}'
        f'{signstr(coef)} ({coef}){rand.choice(consts)})')
    return (deg, ret)

def poly(deg, vars, consts, imag):
    ret = ""
    while i < rand.randrange(deg):
        (degtmp, rtmp) = term(deg, vars, consts, imag)
        i += degtmp
        ret += rtmp

def rational(deg, vars, consts, imag):
    return f'{poly(deg, vars, consts, imag)}/({poly(deg, vars, consts, imag)})'

print(rational(degree, vars, consts, imag))
