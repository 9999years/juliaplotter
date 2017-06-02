# sign funcs
import math
# cli arguments
import argparse
# random rationals
import random

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
    description='Generate a random rational function.'
)

parser.add_argument('--degree', '-d', metavar='degree', type=int,
    default=5, help='The polynomial’s maximum degree.')

parser.add_argument('--variables', '-v', metavar='var', type=str,
    nargs='+', default=['x'], help='An array of variables to use.')

parser.add_argument('--constants', '-c', metavar='const', type=str,
    nargs='+', default=['c'], help='An array of constants to use.')

parser.add_argument('--imaginary', '-i', action='store_true',
    help='Use imaginary numbers.')

# parse arguments, extract variables
args = parser.parse_args()

degree    = args.degree
variables = args.variables
consts    = args.constants
imag      = args.imaginary

consts += ['']

def get_coef(imag):
    max_coef = 20
    coef = (random.random() - 0.5) * max_coef
    if imag:
        coef += get_coef(False) * 1j
    return coef

def term(deg, variables, consts, imag):
    deg = random.randint(0, deg)
    ret = '('
    var = random.choice(variables)
    coef = get_coef(imag)
    if deg != 0:
        ret += f'{var}'
        if deg > 1:
            ret += f'^{deg} '
        else:
            ret += ' '
        ret += f'{signstr(coef)} {abs(coef):.3f}'
    else:
        ret += f'{coef:.3f}'
    ret += f'{random.choice(consts)})'
    return (deg, ret)

def poly(deg, variables, consts, imag):
    ret = ""
    i = 0
    while i < random.randrange(1, deg):
        (degtmp, rtmp) = term(deg, variables, consts, imag)
        i += degtmp
        ret += rtmp
    return ret

def rational(deg, variables, consts, imag):
    ret = poly(deg, variables, consts, imag)
    if random.random() > 0.3:
        ret += f'/({poly(deg, variables, consts, imag)})'
    return ret


print(rational(degree, variables, consts, imag))
